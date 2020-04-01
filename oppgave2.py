# -*- coding: utf-8 -*-

__author__ = "Kjersti Rustad Kvisberg & Ida Lunde Naalsund"
__email__ = "kjkv@nmbu.no, idna@nmbu.no"

import numpy
import pandas as pd
import matplotlib.pyplot as plt

# Defining global variables
h_bar = 1.055E-34
m = 9.109E-31


def psi(pos, pos_0, E, sigma):
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
        return prev_phi + (1j*h_bar/2*m) * \
               (psi(pos + delta_pos, pos_0, E, sigma) -
                2 * psi(pos, pos_0, E, sigma) +
                psi(pos - delta_pos, pos_0, E, sigma)) / delta_pos ** 2 \
               + (V(pos, x_start_step, V_0) *
                  psi(pos, pos_0, E, sigma)) / (1j * h_bar)
    else:
        raise ValueError("Invalid position, must be between 0 and L.")


if __name__ == '__main__':
    sigma = 1E-8
    E = 0.2 * 1.6602E-13
    V_0 = 0.16 * 1.6602E-13
    x_start_step = 100E-9
    x_0 = 50E-9
    L = 200E-9

    delta_pos = 1.5E-10
    delta_t = 2.25E-19
    end_time = 1E-18

    x_values = []
    t_values = []
    phi_values = []

    x = x_0
    while x < L:
        x_values.append(x)
        x += delta_pos
    t = 0
    while t < end_time:
        t_values.append(t)
        t += delta_t

    phi_values = numpy.zeros((len(t_values), len(x_values)))
    prev_phi = 0

    for n in range(0, len(x_values)):
        for m in range(0, len(t_values)):
            x = x_values[n]
            current_phi = phi(
                x, x_0, delta_pos, E, sigma, x_start_step, V_0, prev_phi)
            prev_phi = current_phi
            phi_values.itemset((m, n), current_phi)

    phi_value_data = pd.DataFrame(phi_values, columns=x_values)
    V_values = [V(x, x_start_step, V_0) for x in x_values]

    for i in range(len(t_values)):
        plt.plot(phi_value_data.iloc[i])
        plt.plot(x_values, V_values)
        plt.show()
