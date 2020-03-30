# -*- coding: utf-8 -*-

__author__ = "Kjersti Rustad Kvisberg & Ida Lunde Naalsund"
__email__ = "kjkv@nmbu.no, idna@nmbu.no"

import numpy

# Defining global variables
h_bar = 1.055E-34
m = 9.109E-31


def phi_value(pos, pos_0, E, sigma):
    """ Time-independent wave function for wave packet. """
    k = numpy.sqrt(2*m*E)/h_bar
    a = 1/(numpy.pi**(1/4)*numpy.sqrt(sigma))
    b = (-1/2) * ((pos - pos_0)**2)/(sigma**2)
    c = 1j * k * (pos - pos_0)
    return a * numpy.exp(b) * numpy.exp(c)


def V(x, x_start_step, V_0):
    """Returns potential for given x value."""
    if x > x_start_step:
        return V_0
    else:
        return 0


def phi_zero(pos, pos_0, E, sigma, x_start_step, V_0):
    """ First phi value    """
    delta_pos = 1.5E-10
    return ((1j * h_bar)/(2*m*delta_pos**2) *
            (phi_value(pos+delta_pos, pos_0, E, sigma)
             - 2 * phi_value(pos, pos_0, E, sigma) +
             phi_value(pos - delta_pos, pos_0, E, sigma))
            + (V(pos, x_start_step, V_0) + 1) *
            phi_value(pos, pos_0, E, sigma))


if __name__ == '__main__':
    sigma = 1E-8
    E = 0.2 * 1.6602E-13
    V_0 = 0.16 * 1.6602E-13
    x_start_step = 100E-9
    x_0 = 50E-9
    L = 200E-9
    print(V(125E-9, x_start_step, V_0))
    print(phi_zero(100E-9, x_0, E, sigma, x_start_step, V_0))
