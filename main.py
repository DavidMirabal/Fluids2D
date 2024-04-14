import os
import sys
from tqdm import tqdm
import configparser

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

import utils
import ci
import scheme
import dmb_figures as figures

if len(sys.argv) != 2:
    print("Uso: python main.py config_file")
    sys.exit(1)

config_file = str(sys.argv[1])

# SE LEE EL ARCHIVO DE CONFIGURACIÓN
config = configparser.ConfigParser()
with open(config_file, encoding='utf-8') as file:
    config.read_file(file)

# DIRECTORIOS
nombre = config['DIRECTORIOS']['nombre']
dir_images = config['DIRECTORIOS']['dir_images']
dir_stat = config['DIRECTORIOS']['dir_stat']

if not (os.path.exists(dir_images)):
    os.makedirs(dir_images)
if not (os.path.exists(dir_stat)):
    os.makedirs(dir_stat)

# VISUAL
iteraciones = config.getint('VISUAL', 'iteraciones')
t_show = config.getfloat('VISUAL', 't_show')
dt_show = config.getfloat('VISUAL', 'dt_show')
vmin_rho, vmax_rho = config.getfloat('VISUAL', 'vmin_rho'), config.getfloat(
    'VISUAL', 'vmax_rho'
)
vmin_v, vmax_v = config.getfloat('VISUAL', 'vmin_v'), config.getfloat(
    'VISUAL', 'vmax_v'
)
vmin_p, vmax_p = config.getfloat('VISUAL', 'vmin_p'), config.getfloat(
    'VISUAL', 'vmax_p'
)
d_v = config.getint('VISUAL', 'd_v')
next_show = dt_show + t_show

# SIMULACIÓN
Nx, Ny = config.getint('SIMULACION', 'Nx'), config.getint('SIMULACION', 'Ny')
delta_x = 1 / (Nx - 2)
delta_y = 1 / (Ny - 2)
x_0, x_f = config.getfloat('SIMULACION', 'x_0'), config.getfloat('SIMULACION', 'x_f')
y_0, y_f = config.getfloat('SIMULACION', 'y_0'), config.getfloat('SIMULACION', 'y_f')
rho_0 = config.getfloat('SIMULACION', 'rho_0')
gamma = config.getfloat('SIMULACION', 'gamma')

# CONDICIONES INICIALES
if config['CI']['tipo_ci'] == 'SoundWaves':
    A_rho = config.getfloat('CI', 'A_rho')
    kx, ky = config.getfloat('CI', 'kx'), config.getfloat('CI', 'ky')
    modo_x, modo_y = config.getfloat('CI', 'modo_x'), config.getfloat('CI', 'modo_y')
    P_0 = config.getfloat('CI', 'P_0')
    phi = config.getfloat('CI', 'phi')

    cond_ini = ci.soundwaves_iso(
        (x_0, x_f), (y_0, y_f), delta_x, delta_y, (kx, ky), phi, gamma
    )
    r = cond_ini.r
    X, Y = cond_ini.X, cond_ini.Y
    rho = cond_ini.rho(rho_0, A_rho)
    P = cond_ini.P(P_0)
    c_s = cond_ini.c_s()
    vx, vy, v2, v = cond_ini.v((modo_x, modo_y))
    flux_vx, flux_vy = cond_ini.flux_v()
    rho_e = cond_ini.rho_e()

elif config['CI']['tipo_ci'] == 'GreshoVortex':
    cond_ini = ci.gresho_vortex((x_0, x_f), (y_0, y_f), delta_x, delta_y)
    r = cond_ini.r
    X, Y = cond_ini.X, cond_ini.Y
    rho = cond_ini.rho(rho_0)
    vx, vy, v2, v = cond_ini.v()
    flux_vx, flux_vy = cond_ini.flux_v()
    P = cond_ini.P()
    rho_e = cond_ini.rho_e(gamma)

elif config['CI']['tipo_ci'] == 'GreshoVortexMOD':
    cond_ini = ci.gresho_vortexMOD((x_0, x_f), (y_0, y_f), delta_x, delta_y)
    r = cond_ini.r
    X, Y = cond_ini.X, cond_ini.Y
    rho = cond_ini.rho(rho_0)
    vx, vy, v2, v = cond_ini.v()
    flux_vx, flux_vy = cond_ini.flux_v()
    P = cond_ini.P()
    rho_e = cond_ini.rho_e(gamma)

else:
    sys.exit(
        "Error: El tipo de condición inicial ('tipo_ci') especificado en el archivo de configuración no está disponible."
    )

# TRAZADORES
N_part = config.get('TRAZADORES', 'N_part')
N_part = [int(num) for num in N_part.split(',')]
radios = config.get('TRAZADORES', 'radios')
radios = [float(num) for num in radios.split(',')]
trazadores = [utils.trac(n, r) for n, r in zip(N_part, radios)]

# FIGURAS (condiciones iniciales)
fig_cl = figures.Figura(ratio=2.7)
fig, ax = fig_cl.axs(ncols=3)
im_rho = ax[0].imshow(
    rho,
    extent=(x_0 - delta_x, x_f + delta_x, x_0 - delta_y, x_f + delta_y),
    vmin=vmin_rho,
    vmax=vmax_rho,
)
for k in range(len(trazadores)):
    ax[0].scatter(trazadores[k][:, 0], trazadores[k][:, 1], c='black', s=1)

im_v = ax[1].imshow(
    v,
    extent=(x_0 - delta_x, x_f + delta_x, x_0 - delta_y, x_f + delta_y),
    vmin=vmin_v,
    vmax=vmax_v,
)

ax[1].quiver(
    X[::d_v, ::d_v],
    Y[::d_v, ::d_v],
    vx[::d_v, ::d_v],
    vy[::d_v, ::d_v],
    color='black',
)

im_p = ax[2].imshow(
    P,
    extent=(x_0 - delta_x, x_f + delta_x, x_0 - delta_y, x_f + delta_y),
    vmin=vmin_p,
    vmax=vmax_p,
)

cbar_rho = fig.colorbar(
    im_rho, ax=ax[0], orientation='horizontal', location='top', fraction=0.1
)
cbar_rho.formatter.set_powerlimits(
    (-3, 3),
)
cbar_rho.ax.tick_params(labelsize=10)
cbar_rho.update_ticks()


cbar_v = fig.colorbar(
    im_v, ax=ax[1], orientation='horizontal', location='top', fraction=0.1
)
cbar_v.formatter.set_powerlimits((-3, 3))
cbar_v.ax.tick_params(labelsize=10)
cbar_v.update_ticks()

cbar_p = plt.colorbar(
    im_p, ax=ax[2], orientation='horizontal', location='top', fraction=0.1
)
cbar_p.formatter.set_powerlimits((-3, 3))
cbar_p.ax.tick_params(labelsize=10)
cbar_p.update_ticks()

cbar_rho.set_label(r'$\rho$')
cbar_v.set_label(r'${\rm |\vec{v}|}$')
cbar_p.set_label(r'$P$')

for axis in ax:
    axis.set_xlabel(r'$X$')
    axis.set_ylabel(r'$Y$')
    axis.set_xlim(x_0, x_f)
    axis.set_ylim(y_0, y_f)
    axis.yaxis.set_minor_locator(AutoMinorLocator())
    axis.xaxis.set_minor_locator(AutoMinorLocator())
fig.tight_layout()

ax[1].text(
    0.5,
    -0.45,
    f'$t= {t_show:.2f}$',
    ha='center',
    va='center',
    fontsize=12,
    transform=ax[1].transAxes,
)
fig.savefig(f'{dir_images}im_0.png', dpi=500)

# Momentos y energía
L_angular = np.empty(iteraciones + 1)
P_linear_x = np.empty_like(L_angular)
P_linear_y = np.empty_like(L_angular)
E_total = np.empty_like(L_angular)
time_show_list = np.empty_like(L_angular)

theta_angular = np.arctan2(v, r)
L_angular[0] = np.sum(r * v * rho * np.sin(theta_angular))
P_linear_x[0] = np.sum(rho * vx)
P_linear_y[0] = np.sum(rho * vy)
E_total[0] = np.sum(rho_e)
time_show_list[0] = t_show


def calculate(i):
    global t_show, next_show
    global rho, flux_vx, flux_vy, rho_e
    global P, v, v2, vx, vy, trazadores
    global im_rho, im_v, im_p
    global L_angular, P_linear, E_total, time_show_list

    for axis in ax:
        axis.clear()
        axis.set_xlabel(r'$X$')
        axis.set_ylabel(r'$Y$')
        axis.set_xlim(x_0, x_f)
        axis.set_ylim(y_0, y_f)
        axis.yaxis.set_minor_locator(AutoMinorLocator())
        axis.xaxis.set_minor_locator(AutoMinorLocator())

    while t_show < next_show:
        c_s = np.sqrt(gamma * P / rho)
        delta_t = 0.9 * delta_x / np.max(np.abs(c_s + v))

        LF_actual = scheme.LaxFried2D(
            rho, flux_vx, flux_vy, vx, vy, rho_e, P, gamma, delta_x, delta_y, delta_t
        )
        LF_actual.bound_cond()

        rho, flux_vx, flux_vy, rho_e, vx, vy, v, v2, P = LF_actual.values()

        for traza in trazadores:
            interp_trazadores = utils.semi_interp(X, Y, vx, vy, traza)
            traza[:, 0] = traza[:, 0] + interp_trazadores[0] * delta_t
            traza[:, 1] = traza[:, 1] + interp_trazadores[1] * delta_t

        t_show = t_show + delta_t

    next_show = t_show + dt_show

    im_rho = ax[0].imshow(
        rho,
        extent=(x_0 - delta_x, x_f + delta_x, x_0 - delta_y, x_f + delta_y),
        vmin=vmin_rho,
        vmax=vmax_rho,
    )
    for k in range(len(trazadores)):
        ax[0].scatter(trazadores[k][:, 0], trazadores[k][:, 1], c='black', s=1)

    im_v = ax[1].imshow(
        v,
        extent=(x_0 - delta_x, x_f + delta_x, x_0 - delta_y, x_f + delta_y),
        vmin=vmin_v,
        vmax=vmax_v,
    )

    ax[1].quiver(
        X[::d_v, ::d_v],
        Y[::d_v, ::d_v],
        vx[::d_v, ::d_v],
        vy[::d_v, ::d_v],
        color='black',
    )

    im_p = ax[2].imshow(
        P,
        extent=(x_0 - delta_x, x_f + delta_x, x_0 - delta_y, x_f + delta_y),
        vmin=vmin_p,
        vmax=vmax_p,
    )

    cbar_rho.update_normal(im_rho)
    cbar_v.update_normal(im_v)
    cbar_p.update_normal(im_p)

    theta_angular = np.arctan2(v, r)
    L_angular[i + 1] = np.sum(r * v * rho * np.sin(theta_angular))
    P_linear_x[i + 1] = np.sum(rho * vx)
    P_linear_y[i + 1] = np.sum(rho * vy)
    E_total[i + 1] = np.sum(rho_e)
    time_show_list[i + 1] = next_show

    ax[1].text(
        0.5,
        -0.45,
        f'$t= {t_show:.2f}$',
        ha='center',
        va='center',
        fontsize=12,
        transform=ax[1].transAxes,
    )

    fig.savefig(f'{dir_images}im_{i+1}.png', dpi=500)

    plt.close()


for i in tqdm(range(iteraciones)):
    calculate(i)


fig_cl = figures.Figura()
fig, ax = fig_cl.axs()
ax[0].plot(time_show_list, L_angular, '-', c='black')
ax[0].set_xlabel(r'${\rm Tiempo}$', fontsize=18)
ax[0].set_ylabel(r'${\rm Momento\quad angular}$', fontsize=18)
ax[0].yaxis.set_minor_locator(AutoMinorLocator())
fig.tight_layout()
fig.savefig(f'{dir_stat}{nombre}_mom_angular.pdf')

fig_cl = figures.Figura()
fig, ax = fig_cl.axs()
ax[0].plot(time_show_list, P_linear_x, '-', c='blue', label=r'$p_x$')
ax[0].plot(time_show_list, P_linear_y, '-', c='red', label=r'$p_y$')
ax[0].plot(
    time_show_list,
    np.sqrt(P_linear_x**2 + P_linear_y**2),
    '-',
    c='black',
    label=r'$|\vec{p}|$',
)

ax[0].set_xlabel(r'${\rm Tiempo}$', fontsize=18)
ax[0].set_ylabel(r'${\rm Momento\quad Lineal}$', fontsize=18)
leg = ax[0].legend(fancybox=False, fontsize=14)
leg.get_frame().set_edgecolor('black')
ax[0].yaxis.set_minor_locator(AutoMinorLocator())
fig.tight_layout()
fig.savefig(f'{dir_stat}{nombre}_mom_lineal.pdf')

fig_cl = figures.Figura()
fig, ax = fig_cl.axs()
ax[0].plot(time_show_list, E_total, '-', c='black')
# ax[0].set_ylim(8.9e5, 9e5)
ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
ax[0].set_xlabel(r'${\rm Tiempo}$', fontsize=18)
ax[0].set_ylabel(r'${\rm Energ\acute{i}a\quad Lineal}$', fontsize=18)
ax[0].yaxis.set_minor_locator(AutoMinorLocator())
fig.tight_layout()
fig.savefig(f'{dir_stat}{nombre}_e_total.pdf')


np.savetxt(
    f'{dir_stat}{nombre}_mom_angular.csv', (time_show_list, E_total), delimiter=','
)
