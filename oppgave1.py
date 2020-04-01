# -*- coding: utf-8 -*-

__author__ = "Kjersti Rustad Kvisberg & Ida Lunde Naalsund"
__email__ = "kjkv@nmbu.no, idna@nmbu.no"

# Problem 1

import numpy as np
import matplotlib.pyplot as plt

# Defining global variables
h_bar = 1.055E-34
m = 9.109E-31


def k(E, m):
    """ Calculates wave number k. """
    return np.sqrt(2*m*E)/h_bar


def k_0(E, V_0, m):
    """ Calculated wave number k_0 for step potential. """
    return np.sqrt(2*m*(E-V_0))/h_bar


def trans_coef(E, V_0):
    """ Calculates transmission coefficient for particle interacting with
    step potential. """
    return (4*k(E, m)*k_0(E, V_0, m))/(k(E, m) + k_0(E, V_0, m))**2


def ref_coef(E, V_0):
    """ Calculates reflection coefficient for particle interacting with step
    potential. """
    return (k(E, m) - k_0(E, V_0, m))**2 / (k(E, m) + k_0(E, V_0, m))**2


if __name__ == '__main__':
    E = 0.20 * 1.6602E-13
    V_0 = 0.16 * 1.6602E-13

    T = trans_coef(E, V_0)
    R = ref_coef(E, V_0)
    print(f' The transmission coefficient is {T:.3f}')
    print(f' The reflection coefficient is {R:.3f}')
    print(f' The sum of the transmission and reflection coefficients is '
          f'{(R + T):.3f}')

    E = np.linspace(0.2*1.6602E-13, 0.4*1.6602E-13, 1000)
    plt.plot(E, trans_coef(E, V_0))
    plt.title('Transmission coefficient for different particle energies')
    plt.xlabel('Energy [J]')
    plt.ylabel('Transmission coefficient')
    plt.show()


