# -*- coding: utf-8 -*-

__author__ = "Ida Lunde Naalsund & Kjersti Rustad Kvisberg"
__email__ = "idna@nmbu.no, kjkv@nmbu.no"

import numpy as np
import matplotlib.pyplot as plt

# Defining global variables
h_bar = 1.055E-34
m = 9.109E-31

def k(E):
    """ Calculates wave number k. """
    return np.sqrt(2*m*E)/h_bar


def k_0(E, V_0):
    """ Calculates wave number k_0 for step potential. """
    return np.sqrt(2*m*(E-V_0))/h_bar

def k_x(E, alpha):
    """ Calculates x component of wave number k"""
    return k(E)*np.cos(alpha)

def k_y(E, alpha):
    """ Calculates y component of wave number k"""
    return k(E)*np.sin(alpha)

def k_0x(E, V_0, alpha):
    """ Calculates x component of wave number k_0"""
    return np.sqrt(k_0(E, V_0)**2 - k_y(E, alpha)**2)


def transmission_coef(E, V_0, alpha):
    return np.sqrt((4*k_x(E, alpha)*k_0x(E, V_0, alpha))/
                   (k_x(E, alpha) + k_0x(E, V_0, alpha))**2)**2


if __name__ == "__main__":
    E = np.linspace(0.2, 0.4, 100)*1.6022E-19
    alpha = np.linspace(0, np.pi/9, 100)
    V_0 = 0.16*1.6022E-19

    trans = np.empty((100, 100))
    for energy in E:
        T = transmission_coef(energy, V_0, alpha)
        T.reshape((1, 100))
        trans = np.append(trans, np.array(T), axis=0)

    print(trans)


