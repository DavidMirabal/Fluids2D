import numpy as np
import utils


def LFS2D(
    u: np.ndarray, der_flux_u_x: np.ndarray, der_flux_u_y: np.ndarray, delta_t: float
) -> np.ndarray:
    """
    Definición del esquema LaxFriedrichs en 2D para 1 variable conservativa.add()

    :param u: Variable
    :param der_flux_u_x: Derivada de la variable en la dirección X
    :param der_flux_u_y: Derivada de la variable en la dirección Y
    :param delta_t: Espaciado temporal
    """
    return 0.25 * (
        u[1:-1, 2:] + u[1:-1, :-2] + u[2:, 1:-1] + u[:-2, 1:-1]
    ) - 0.5 * delta_t * (der_flux_u_x + der_flux_u_y)


class LaxFried2D:
    def __init__(
        self,
        rho: np.ndarray,
        flux_vx: np.ndarray,
        flux_vy: np.ndarray,
        vx: np.ndarray,
        vy: np.ndarray,
        rho_e: np.ndarray,
        P: np.ndarray,
        gamma: float,
        delta_x: float,
        delta_y: float,
        delta_t: float,
    ):
        """
        Calcula el siguiente paso temporal de rho, flux_v, rho_e con el esquema LaxFriedrichs en 2D.

        :param rho: Densidad de materia
        :param flux_vx: Flujo de materia en la dirección X
        :param flux_vy: Flujo de materia en la dirección Y
        :param vx: Velocidad en la dirección X
        :param vy: Velocidad en la dirección Y
        :param rho_e: Densidad de energía
        :param P: Presión
        :param gamma: Índice adiabático
        :param delta_x: Espaciado en X del grid
        :param delta_y: Espaciado en Y del grid
        :param delta_t: Espaciado en tiempo (paso de integración)
        """

        self.rho = rho
        self.flux_vx, self.flux_vy = flux_vx, flux_vy
        self.vx, self.vy = vx, vy
        self.rho_e = rho_e
        self.P = P
        self.gamma = gamma

        der_flux_x = utils.der(flux_vx).euler_medio2d(delta_x, dim='i')
        der_flux_y = utils.der(flux_vy).euler_medio2d(delta_y, dim='j')

        der_mom_vx_x = utils.der(flux_vx * vx + P).euler_medio2d(delta_x, dim='i')
        der_mom_vcross = utils.der(flux_vx * vy).euler_medio2d(
            dim='both', delta_x=delta_x, delta_y=delta_y
        )
        der_mom_vx_y = der_mom_vcross[1]
        der_mom_vy_x = der_mom_vcross[0]
        der_mom_vy_y = utils.der(flux_vy * vy + P).euler_medio2d(delta_y, dim='j')

        rho_e_P = rho_e + P
        der_ene_x = utils.der(rho_e_P * vx).euler_medio2d(delta_x, dim='i')
        der_ene_y = utils.der(rho_e_P * vy).euler_medio2d(delta_y, dim='j')

        rho_n1 = LFS2D(rho, der_flux_x, der_flux_y, delta_t)

        flux_vx_n1 = LFS2D(flux_vx, der_mom_vx_x, der_mom_vx_y, delta_t)
        flux_vy_n1 = LFS2D(flux_vy, der_mom_vy_x, der_mom_vy_y, delta_t)

        rho_e_n1 = LFS2D(rho_e, der_ene_x, der_ene_y, delta_t)

        self.rho[1:-1, 1:-1] = rho_n1
        self.flux_vx[1:-1, 1:-1] = flux_vx_n1
        self.flux_vy[1:-1, 1:-1] = flux_vy_n1
        self.rho_e[1:-1, 1:-1] = rho_e_n1

    def bound_cond(self):
        """
        Condiciones de contorno periódicas.
        """
        for var in [self.rho, self.flux_vx, self.flux_vy, self.rho_e]:
            var[0, :] = var[-2, :]
            var[-1, :] = var[1, :]

            var[:, 0] = var[:, -2]
            var[:, -1] = var[:, 1]

    def values(self) -> tuple:
        """
        Devuelve las variables calculadas

        :return: Densidad de materia, flujos de materia, densidad de energía, velocidades y presión
        """
        vx = self.flux_vx / self.rho
        vy = self.flux_vy / self.rho
        v2 = vx * vx + vy * vy
        v = np.sqrt(v2)
        P = (self.rho_e - 0.5 * self.rho * v2) * (self.gamma - 1)

        return self.rho, self.flux_vx, self.flux_vy, self.rho_e, vx, vy, v, v2, P
