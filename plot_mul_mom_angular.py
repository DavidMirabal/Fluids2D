import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

import utils
import dmb_figures as figures

# Directorio con los csv generados por main.py
dir_stats = 'vortex/stats/'
archivos_stats = [f for f in os.listdir(dir_stats) if f.endswith('csv')]


try:
    archivos_stats = sorted(archivos_stats, key=utils.extraer_num)
except:
    pass

times = [None] * len(archivos_stats)
moms_angular = [None] * len(times)
dif_mom_angular = [None] * len(times)
dif_times = np.empty(len(archivos_stats))
for i, archiv in enumerate(archivos_stats):
    simulation = np.loadtxt(f'{dir_stats}{archiv}', delimiter=',')
    times[i] = simulation[0]
    moms_angular[i] = simulation[1]
    dif_times[i] = times[i][-1] - times[i][0]
    dif_mom_angular[i] = np.abs(1 - moms_angular[i] / np.max(moms_angular[i]))

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage[varg]{txfonts}'
plt.rcParams['text.antialiased'] = True
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13

fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

colors = figures.generar_colores(len(times))

nombres = [os.path.splitext(f)[0].split('_mom_angular')[0] for f in archivos_stats]
for i in range(len(times)):
    ax[0].plot(
        times[i],
        moms_angular[i] / np.max(moms_angular[i]),
        '-',
        c=colors[i],
        label=nombres[i],
    )
    ax[1].plot(
        times[i], 1 - moms_angular[i] / np.max(moms_angular[i]), '-', c=colors[i]
    )

ax[1].set_xlabel(r'${\rm Tiempo}$', fontsize=18)
ax[0].set_ylabel(r'${\rm Momento\quad angular\quad (norm.)}$', fontsize=18)
ax[1].set_ylabel(r'$\Delta$', fontsize=18)
ax[0].yaxis.set_minor_locator(AutoMinorLocator())
ax[1].yaxis.set_minor_locator(AutoMinorLocator())
leg = ax[0].legend(fancybox=False, fontsize=14)
leg.get_frame().set_edgecolor('black')
fig.tight_layout()
fig.subplots_adjust(hspace=0)
fig.savefig(f'{dir_stats}MUL_mom_angular.pdf')


try:
    fig_cl = figures.Figura()
    fig, ax = fig_cl.axs()

    numeros = [int(nom.split('N')[-1]) for nom in nombres]
    ax[0].plot(
        numeros,
        np.sum(dif_mom_angular, axis=1) ** 2 / dif_times,
        ls='-',
        c='green',
        marker='o',
        ms=5,
    )
    ax[0].set_yscale('log')
    ax[0].set_xlabel(r'${\rm Puntos\quad por\quad dimensi\acute{o}n}$', fontsize=18)
    ax[0].set_ylabel(r'${\rm Par\acute{a}metro\quad de\quad error}$', fontsize=18)
    fig.savefig(f'{dir_stats}CONV_mom_angular.pdf')


except:
    pass
