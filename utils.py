import numpy as np


class der:
    def __init__(self, f):
        self.f = f

    def euler(self, h):
        df = (self.f[1:] - self.f[:-1]) / h
        return df

    def euler_medio(self, h):
        df = (self.f[2:] - self.f[:-2]) / (2 * h)
        return df

    def euler_medio2d(
        self,
        h: float = None,
        dim: str = 'i',
        delta_x: float = None,
        delta_y: float = None,
    ) -> np.ndarray:
        """
        Calcula la derivada euleriana en 2D

        :param h: Espaciado en X/Y
        :param dim: Dimensión sobre la que se deriva (i==X, j==Y). Si se selecciona both se calculan las 2
        :param delta_x: Espaciado en X (solo para dim=both).
        :param delta_y: Espaciado en Y (solo para dim=both)
        """

        if dim == 'i':
            df = (self.f[1:-1, 2:] - self.f[1:-1, :-2]) / (2 * h)
        elif dim == 'j':
            df = (self.f[2:, 1:-1] - self.f[:-2, 1:-1]) / (2 * h)
        elif dim == 'both':
            df = np.array(
                [
                    (self.f[1:-1, 2:] - self.f[1:-1, :-2]) / (2 * delta_x),
                    (self.f[2:, 1:-1] - self.f[:-2, 1:-1]) / (2 * delta_y),
                ]
            )
        return df

    def L1(self, df_true, df):
        N = len(df)
        return 1 / N * np.sum(np.abs(df_true - df))


def trac(N: int, r: float) -> np.ndarray:
    """
    Calcula la posición de trazadores (partículas sobre el fluido) situados equiespaciadamente a un radio

    :param N: Número de partículas
    :param r: Radio donde situar los trazadores
    """
    theta = np.linspace(0, 2 * np.pi, N, endpoint=False)
    return np.column_stack((r * np.cos(theta), r * np.sin(theta)))


def semi_interp(
    X: np.ndarray, Y: np.ndarray, x_p: np.ndarray, y_p: np.ndarray, trazador: np.ndarray
) -> tuple:
    """
    Asigna el valor de x_p y y_p más cercano a la posición del trazador

    :param X: Valores de X en el grid
    :param Y: Valores de Y en el grid
    :param x_p: Mapa del que se seleccionará el valor más cercano al trazador
    :param y_p: Mapa del que se seleccionará el valor más cercano al trazador
    :param trazador: Posiciones de los trazadores
    """
    index_x = np.empty(len(trazador[:, 0]))
    index_y = np.empty_like(index_x)

    for i, tra in enumerate(trazador):
        index_x[i] = np.argmin(np.abs(X[0, :] - tra[0]))
        index_y[i] = np.argmin(np.abs(Y[:, 0] - tra[1]))

    return (
        x_p[np.array(index_y, dtype=int), np.array(index_x, dtype=int)],
        y_p[np.array(index_y, dtype=int), np.array(index_x, dtype=int)],
    )


def extraer_num(nombre):
    return int(nombre.split('N')[-1].split('_')[0])
