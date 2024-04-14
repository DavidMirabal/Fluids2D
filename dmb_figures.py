import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import numpy as np


def generar_colores(N):
    valores = np.linspace(0, 1, N)
    colores = plt.cm.rainbow(valores)
    colores_hex = [to_hex(color) for color in colores]
    return colores_hex


class Figura:
    def __init__(self, ancho=6.4, ratio=1.33, dpi=200):
        # Las siguientes líneas son para homogeneizar las figuras con el estilo de MNRAS
        # Se puede comentar si no funcionan correctamente los paquetes de LaTeX
        plt.rcParams['text.usetex'] = True
        plt.rcParams['text.latex.preamble'] = r'\usepackage[varg]{txfonts}'
        plt.rcParams['text.antialiased'] = True
        plt.rcParams['xtick.labelsize'] = 13
        plt.rcParams['ytick.labelsize'] = 13
        plt.rcParams['image.cmap'] = 'rainbow'

        fig = plt.figure(figsize=(ancho, ancho / ratio), dpi=dpi)

        self.fig = fig
        self.dpi = dpi
        self.ax = []

    def axs(self, ncols=1, nrows=1):
        k = 0
        for i in range(nrows):
            for j in range(ncols):
                self.ax.append(self.fig.add_subplot(nrows, ncols, k + 1))
                self.ax[-1].tick_params(axis='both')
                k = k + 1

        return self.fig, self.ax
