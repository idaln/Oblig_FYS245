# -*- coding: utf-8 -*-

__author__ = "Ida Lunde Naalsund & Kjersti Rustad Kvisberg"
__email__ = "idna@nmbu.no, kjkv@nmbu.no"

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Defining global variables
h_bar = 1.055E-34
m = 9.109E-31


def k(E):
    """ Calculates wave number k. """
    return np.sqrt(2 * m * E) / h_bar


def k_0(E, V_0):
    """ Calculates wave number k_0 for step potential. """
    return np.sqrt(2 * m * (E - V_0)) / h_bar


def k_x(E, alpha):
    """ Calculates x component of wave number k"""
    return k(E) * np.cos(alpha)


def k_y(E, alpha):
    """ Calculates y component of wave number k"""
    return k(E) * np.sin(alpha)


def k_0x(E, V_0, alpha):
    """ Calculates x component of wave number k_0"""
    return np.sqrt(k_0(E, V_0) ** 2 - k_y(E, alpha) ** 2)


def transmission_coef(E, V_0, alpha):
    """ Calculates """
    return np.sqrt(((4 * k_x(E, alpha) * k_0x(E, V_0, alpha)) /
                   (k_x(E, alpha) + k_0x(E, V_0, alpha)) ** 2) ** 2)


if __name__ == "__main__":
    n = 100
    E = np.linspace(0.2, 0.4, n) * 1.6022E-19
    alpha = np.linspace(0, np.pi/9, n)
    V_0 = 0.16 * 1.6022E-19

    trans_coef = np.empty((0, n))
    for energy in E:
        T = transmission_coef(energy, V_0, alpha)
        trans_coef = np.vstack((trans_coef, np.array(T)))

    num_ticks = np.int(n/10)
    xticks = np.linspace(0, n - 1, num_ticks, dtype=np.int)
    yticks = np.linspace(0, n - 1, num_ticks, dtype=np.int)
    xticklabels = [round(alpha[idx] * (180 / np.pi), 1) for idx in xticks]
    yticklabels = [round(E[idx]/1.6022E-19, 3) for idx in yticks]

    sns.set()
    ax = sns.heatmap(
        trans_coef, cmap="GnBu",
        xticklabels=xticklabels,  yticklabels=yticklabels
    )
    ax.set_yticks(yticks)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=0)
    plt.xlabel("Incident angle [degrees]")
    plt.ylabel("Energy values [eV]")
    plt.title("Transmission coefficient of incident angle and energy")
    plt.show()

    plt.plot(alpha*180/np.pi, transmission_coef(E[-1], V_0, alpha))
    plt.xlabel("Alpha values [degrees]")
    plt.ylabel("Transmission coefficient")
    plt.title("Transmission coefficient of incident angle with energy 0.40 eV")
    plt.show()

    print(k_x(E=0.2 * 1.6022E-19, alpha=alpha))
    print(k_0x(E=0.2 * 1.6022E-19, V_0=V_0, alpha=alpha))
    print(k_y(E=0.2 * 1.6022E-19, alpha=alpha))