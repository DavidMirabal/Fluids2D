import numpy as np


# Condiciones iniciales


# Onda de sonido en un medio isotermo
class soundwaves_iso:  # es necesario que ejecutar funciones de esta clase en el orden en el que aparecen
    def __init__(
        self,
        x_0_f: tuple,
        y_0_f: tuple,
        delta_x: float,
        delta_y: float,
        k: tuple,
        phi: float,
        gamma: float,
    ):
        """
        Calcula el grid.

        :param x_0_f: Intervalo del dominio (en X) en el que se calcula
        :param y_0_f: Intervalo del dominio (en Y) en el que se calcula
        :param delta_x: Espaciado en X del grid
        :param delta_y: Espaciado en Y del grid
        :param k: Número de onda para ambos ejes
        :param phi: Fase de la onda
        :param gamma: Índice adiabático
        """

        x = np.arange(x_0_f[0] - delta_x, x_0_f[1] + delta_x, delta_x)
        y = np.arange(y_0_f[0] - delta_y, y_0_f[1] + delta_y, delta_y)

        self.delta_x = delta_x
        self.delta_y = delta_y

        self.X, self.Y = np.meshgrid(x, y)

        self.r2 = self.X * self.X + self.Y * self.Y
        self.r = np.sqrt(self.r2)

        self.k = k
        self.phi = phi
        self.gamma = gamma

    def rho(self, rho_0: float, A_rho: float) -> np.ndarray:
        """
        Crea un mapa de densidad constante en el espacio.

        :param rho_0: Valor promedio del mapa de densidad
        :param A_rho: Valor del pico de la perturbación
        :return: Mapa de densidad
        """
        self.rho_0, self.A_rho = rho_0, A_rho
        self.A_P = self.gamma * self.A_rho
        self.rho = self.rho_0 + self.A_rho * np.cos(
            self.k[0] * self.X + self.k[1] * self.Y + self.phi
        )
        return self.rho

    def P(self, P_0: float) -> np.ndarray:
        """
        Calcula el mapa de presiones.

        :return: Mapa de presiones
        """

        self.P_0 = P_0
        self.P = self.P_0 + self.A_P * np.cos(
            self.k[0] * self.X + self.k[1] * self.Y + self.phi
        )

        return self.P

    def c_s(self) -> float:
        """
        Calcula la velocidad del sonido (sin contar la perturbación).

        :return: Mapa de presiones
        """
        self.c_s = np.sqrt(self.gamma * self.P_0 / self.rho_0)
        return self.c_s

    def v(self, modo: tuple = (+1, +1)) -> tuple:
        """
        Calcula el mapa de velocidades.

        :param modo: Modo positivo o negativo de velocidad en X y en Y
        :return: Mapa de velocidad, modulo de velocidades al cuadrado y modulo de velocidades
        """
        A_v = np.empty(2)
        for i, mod in enumerate(modo):
            A_v[i] = mod * self.c_s * self.A_P

        self.vx = A_v[0] * np.cos(self.k[0] * self.X + self.k[1] * self.Y + self.phi)
        self.vy = A_v[1] * np.cos(self.k[0] * self.X + self.k[1] * self.Y + self.phi)

        self.v2 = self.vx * self.vx + self.vy * self.vy
        v = np.sqrt(self.v2)

        return self.vx, self.vy, self.v2, v

    def flux_v(self) -> tuple:
        """
        Calcula el mapa de flujos de densidad.

        :return: Mapa de flujos de densidad
        """

        return self.rho * self.vx, self.rho * self.vy

    def rho_e(self) -> np.ndarray:
        """
        Calcula el mapa de densidades de energía.

        :return: Mapa de densidades de energía
        """

        return self.P / (self.gamma - 1) + 0.5 * self.rho * self.v2


# Gresho Vortex
class gresho_vortex:  # es necesario que ejecutar funciones de esta clase en el orden en el que aparecen
    def __init__(self, x_0_f: tuple, y_0_f: tuple, delta_x: float, delta_y: float):
        """
        Calcula el grid.

        :param x_0_f: Intervalo del dominio (en X) en el que se calcula
        :param y_0_f: Intervalo del dominio (en Y) en el que se calcula
        :param delta_x: Espaciado en X del grid
        :param delta_y: Espaciado en Y del grid
        """

        x = np.arange(x_0_f[0] - delta_x, x_0_f[1] + delta_x, delta_x)
        y = np.arange(y_0_f[0] - delta_y, y_0_f[1] + delta_y, delta_y)

        self.delta_x = delta_x
        self.delta_y = delta_y

        self.X, self.Y = np.meshgrid(x, y)

        self.theta = np.arctan2(self.Y, self.X)

        self.r2 = self.X * self.X + self.Y * self.Y
        self.r = np.sqrt(self.r2)

        self.r_mask_0p2 = self.r <= 0.2
        self.r_mask_0p2_0p4 = (self.r > 0.2) & (self.r < 0.4)
        self.r_mask_0p4 = self.r >= 0.4

    def rho(self, rho_0: float) -> np.ndarray:
        """
        Crea un mapa de densidad constante en el espacio.

        :param rho_0: Valor del mapa de densidad
        :return: Mapa de densidad
        """

        self.rho = np.ones_like(self.r) * rho_0
        return self.rho

    def v(self) -> tuple:
        """
        Calcula el mapa de velocidades.

        :return: Mapa de velocidad, modulo de velocidades al cuadrado y modulo de velocidades
        """

        vx = np.empty_like(self.r)
        vy = np.empty_like(self.r)

        vx[self.r_mask_0p2] = (
            5 * self.r[self.r_mask_0p2] * np.sin(self.theta[self.r_mask_0p2])
        )
        vy[self.r_mask_0p2] = (
            5 * self.r[self.r_mask_0p2] * -np.cos(self.theta[self.r_mask_0p2])
        )

        vx[self.r_mask_0p2_0p4] = (2 - 5 * self.r[self.r_mask_0p2_0p4]) * np.sin(
            self.theta[self.r_mask_0p2_0p4]
        )
        vy[self.r_mask_0p2_0p4] = (2 - 5 * self.r[self.r_mask_0p2_0p4]) * -np.cos(
            self.theta[self.r_mask_0p2_0p4]
        )

        vx[self.r_mask_0p4] = 0
        vy[self.r_mask_0p4] = 0

        self.vx, self.vy = vx, vy

        self.v2 = vx * vx + vy * vy
        v = np.sqrt(self.v2)

        return self.vx, self.vy, self.v2, v

    def flux_v(self) -> tuple:
        """
        Calcula el mapa de flujos de densidad.

        :return: Mapa de flujos de densidad
        """

        return self.rho * self.vx, self.rho * self.vy

    def P(self) -> np.ndarray:
        """
        Calcula el mapa de presiones.

        :return: Mapa de presiones
        """

        P = np.empty_like(self.r)
        P[self.r_mask_0p2] = 5 + 12.5 * self.r2[self.r_mask_0p2]
        P[self.r_mask_0p2_0p4] = (
            9
            + 12.5 * self.r2[self.r_mask_0p2_0p4]
            - 20 * self.r[self.r_mask_0p2_0p4]
            + 4 * np.log(5 * self.r[self.r_mask_0p2_0p4])
        )
        P[self.r >= 0.4] = 3 + 4 * np.log(2)

        self.P = P

        return self.P

    def rho_e(self, gamma: float) -> np.ndarray:
        """
        Calcula el mapa de densidades de energía.

        :param gamma: Índice adiabático
        :return: Mapa de densidades de energía
        """

        return self.P / (gamma - 1) + 0.5 * self.rho * self.v2


# CI FOR FUN:
# Gresho Vortex MOD (cambio el logaritmo neperiano por logaritmo en base 10)
class gresho_vortexMOD:  # es necesario que ejecutar funciones de esta clase en el orden en el que aparecen
    def __init__(self, x_0_f: tuple, y_0_f: tuple, delta_x: float, delta_y: float):
        """
        Calcula el grid.

        :param x_0_f: Intervalo del dominio (en X) en el que se calcula
        :param y_0_f: Intervalo del dominio (en Y) en el que se calcula
        :param delta_x: Espaciado en X del grid
        :param delta_y: Espaciado en Y del grid
        """

        x = np.arange(x_0_f[0] - delta_x, x_0_f[1] + delta_x, delta_x)
        y = np.arange(y_0_f[0] - delta_y, y_0_f[1] + delta_y, delta_y)

        self.delta_x = delta_x
        self.delta_y = delta_y

        self.X, self.Y = np.meshgrid(x, y)

        self.theta = np.arctan2(self.Y, self.X)

        self.r2 = self.X * self.X + self.Y * self.Y
        self.r = np.sqrt(self.r2)

        self.r_mask_0p2 = self.r <= 0.2
        self.r_mask_0p2_0p4 = (self.r > 0.2) & (self.r < 0.4)
        self.r_mask_0p4 = self.r >= 0.4

    def rho(self, rho_0: float) -> np.ndarray:
        """
        Crea un mapa de densidad constante en el espacio.

        :param rho_0: Valor del mapa de densidad
        :return: Mapa de densidad
        """

        self.rho = np.ones_like(self.r) * rho_0
        return self.rho

    def v(self) -> tuple:
        """
        Calcula el mapa de velocidades.

        :return: Mapa de velocidad, modulo de velocidades al cuadrado y modulo de velocidades
        """

        vx = np.empty_like(self.r)
        vy = np.empty_like(self.r)

        vx[self.r_mask_0p2] = (
            5 * self.r[self.r_mask_0p2] * np.sin(self.theta[self.r_mask_0p2])
        )
        vy[self.r_mask_0p2] = (
            5 * self.r[self.r_mask_0p2] * -np.cos(self.theta[self.r_mask_0p2])
        )

        vx[self.r_mask_0p2_0p4] = (2 - 5 * self.r[self.r_mask_0p2_0p4]) * np.sin(
            self.theta[self.r_mask_0p2_0p4]
        )
        vy[self.r_mask_0p2_0p4] = (2 - 5 * self.r[self.r_mask_0p2_0p4]) * -np.cos(
            self.theta[self.r_mask_0p2_0p4]
        )

        vx[self.r_mask_0p4] = 0
        vy[self.r_mask_0p4] = 0

        self.vx, self.vy = vx, vy

        self.v2 = vx * vx + vy * vy
        v = np.sqrt(self.v2)

        return self.vx, self.vy, self.v2, v

    def flux_v(self) -> tuple:
        """
        Calcula el mapa de flujos de densidad.

        :return: Mapa de flujos de densidad
        """

        return self.rho * self.vx, self.rho * self.vy

    def P(self) -> np.ndarray:
        """
        Calcula el mapa de presiones.

        :return: Mapa de presiones
        """

        P = np.empty_like(self.r)
        P[self.r_mask_0p2] = 5 + 12.5 * self.r2[self.r_mask_0p2]
        P[self.r_mask_0p2_0p4] = (
            9
            + 12.5 * self.r2[self.r_mask_0p2_0p4]
            - 20 * self.r[self.r_mask_0p2_0p4]
            + 4 * np.log10(5 * self.r[self.r_mask_0p2_0p4])
        )
        P[self.r >= 0.4] = 3 + 4 * np.log10(2)

        self.P = P

        return self.P

    def rho_e(self, gamma: float) -> np.ndarray:
        """
        Calcula el mapa de densidades de energía.

        :param gamma: Índice adiabático
        :return: Mapa de densidades de energía
        """

        return self.P / (gamma - 1) + 0.5 * self.rho * self.v2


class Explosion:  # es necesario que ejecutar funciones de esta clase en el orden en el que aparecen
    def __init__(self, x_0_f: tuple, y_0_f: tuple, delta_x: float, delta_y: float):
        """
        Calcula el grid.

        :param x_0_f: Intervalo del dominio (en X) en el que se calcula
        :param y_0_f: Intervalo del dominio (en Y) en el que se calcula
        :param delta_x: Espaciado en X del grid
        :param delta_y: Espaciado en Y del grid
        """

        x = np.arange(x_0_f[0] - delta_x, x_0_f[1] + delta_x, delta_x)
        y = np.arange(y_0_f[0] - delta_y, y_0_f[1] + delta_y, delta_y)

        self.delta_x = delta_x
        self.delta_y = delta_y

        self.X, self.Y = np.meshgrid(x, y)

        self.theta = np.arctan2(self.Y, self.X)

        self.r2 = self.X * self.X + self.Y * self.Y
        self.r = np.sqrt(self.r2)

        self.r_mask_0p2 = self.r <= 0.2
        self.r_mask_0p2_0p4 = (self.r > 0.2) & (self.r < 0.4)
        self.r_mask_0p4 = self.r >= 0.4

    def rho(self, rho_0: float) -> np.ndarray:
        """
        Crea un mapa de densidad constante en el espacio.

        :param rho_0: Valor del mapa de densidad
        :return: Mapa de densidad
        """

        self.rho = np.ones_like(self.r) * rho_0
        return self.rho

    def v(self) -> tuple:
        """
        Calcula el mapa de velocidades.

        :return: Mapa de velocidad, modulo de velocidades al cuadrado y modulo de velocidades
        """

        vx = np.empty_like(self.r)
        vy = np.empty_like(self.r)

        vx[self.r_mask_0p2] = np.cos(self.theta[self.r_mask_0p2])
        vy[self.r_mask_0p2] = np.sin(self.theta[self.r_mask_0p2])

        vx[self.r_mask_0p2_0p4] = (2 - 5 * self.r[self.r_mask_0p2_0p4]) * np.sin(
            self.theta[self.r_mask_0p2_0p4]
        )
        vy[self.r_mask_0p2_0p4] = (2 - 5 * self.r[self.r_mask_0p2_0p4]) * -np.cos(
            self.theta[self.r_mask_0p2_0p4]
        )

        vx[self.r_mask_0p4] = 0
        vy[self.r_mask_0p4] = 0

        self.vx, self.vy = vx, vy

        self.v2 = vx * vx + vy * vy
        v = np.sqrt(self.v2)

        return self.vx, self.vy, self.v2, v

    def flux_v(self) -> tuple:
        """
        Calcula el mapa de flujos de densidad.

        :return: Mapa de flujos de densidad
        """

        return self.rho * self.vx, self.rho * self.vy

    def P(self) -> np.ndarray:
        """
        Calcula el mapa de presiones.

        :return: Mapa de presiones
        """

        P = np.empty_like(self.r)
        P[self.r_mask_0p2] = 5 + 12.5 * self.r2[self.r_mask_0p2]
        P[self.r_mask_0p2_0p4] = (
            9
            + 12.5 * self.r2[self.r_mask_0p2_0p4]
            - 20 * self.r[self.r_mask_0p2_0p4]
            + 4 * np.log10(5 * self.r[self.r_mask_0p2_0p4])
        )
        P[self.r >= 0.4] = 3 + 4 * np.log10(2)

        self.P = P

        return self.P

    def rho_e(self, gamma: float) -> np.ndarray:
        """
        Calcula el mapa de densidades de energía.

        :param gamma: Índice adiabático
        :return: Mapa de densidades de energía
        """

        return self.P / (gamma - 1) + 0.5 * self.rho * self.v2
