# -*- coding: utf-8 -*-

__author__ = "Kjersti Rustad Kvisberg & Ida Lunde Naalsund"
__email__ = "kjkv@nmbu.no, idna@nmbu.no"

import numpy

# Defining global variables
h_bar = 1.055E-34
m = 9.109E-31


def k(E, m):
    """ Calculates wave number k. """
    return np.sqrt(2*m*E)/h_bar


def phi_value(pos, pos_0, energy, sigma):
    """ Time independent wave function for wave packet. """
    k = k(energy, m)
    a = 1/(numpy.pi**(1/4)*numpy.sqrt(sigma))
    b = (-1/2) * ((pos - pos_0)**2)/(sigma**2)
    c = 1j * k * (pos - pos_0)
    return a * numpy.exp(b) * numpy.exp(c)

if __name__ == '__main__':
    