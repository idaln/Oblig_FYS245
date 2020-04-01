# -*- coding: utf-8 -*-

__author__ = "Kjersti Rustad Kvisberg & Ida Lunde Naalsund"
__email__ = "kjkv@nmbu.no, idna@nmbu.no"

import numpy

# Defining global variables
h_bar = 1.055E-34
m = 9.109E-31


def psi_value(pos, pos_0, E, sigma):
    """ Time-independent wave function for wave packet. """
    k = numpy.sqrt(2*m*E)/h_bar
    a = 1/(numpy.pi**(1/4)*numpy.sqrt(sigma))
    b = (-1/2) * ((pos - pos_0)**2)/(sigma**2)
    c = 1j * k * (pos - pos_0)
    if pos == 0 or pos == L:
        return 0
    else:
        return a * numpy.exp(b) * numpy.exp(c)


def V(x, x_start_step, V_0):
    """Returns potential for given x value."""
    if x > x_start_step:
        return V_0
    else:
        return 0


def phi(pos, pos_0, delta_pos, E, sigma, x_start_step, V_0, prev_phi):
    """ Calculates phi value """
    if pos == 0 or pos == L:
        return 0
    elif 0 < pos < L:
        return prev_phi + (1j*h_bar/2*m)*\
            (psi_value(pos+delta_pos, pos_0, E, sigma) -
                2*psi_value(pos, pos_0, E, sigma) +
                psi_value(pos-delta_pos, pos_0, E, sigma))/(delta_pos)**2 + \
            (V(pos, x_start_step, V_0)*
             psi_value(pos, pos_0, E, sigma))/1j*h_bar
    else:
        raise ValueError("Invalid position, must be between 0 and L.")




if __name__ == '__main__':
    sigma = 1E-8
    E = 0.2 * 1.6602E-13
    V_0 = 0.16 * 1.6602E-13
    x_start_step = 100E-9
    x_0 = 50E-9
    L = 200E-9
    #print(psi_value(L, x_0, E, sigma))
    #print(V(125E-9, x_start_step, V_0))
    #print(phi(100E-9, x_0, E, sigma, x_start_step, V_0, 0.00000000003))
    prev_phi = 0
    delta_pos = 1.5E-10
    delta_t = 2.25E-19

    end_time = 5E-18
    x = 0
    x_values = []
    phi_values = []
    while x < L:
        t = 0
        while t < end_time:
            current_phi = (phi(x, x_0, delta_pos, E, sigma, x_start_step,
                              V_0, prev_phi))
            prev_phi = current_phi
            phi_values.append(current_phi)

            t += delta_t

        x += delta_pos




