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
    k = (numpy.sqrt(2*m*E))/h_bar
    a = 1/(numpy.pi**(1/4)*numpy.sqrt(sigma))
    b = (-1/2) * ((pos - pos_0)**2)/(sigma**2)
    c = 1j * k * (pos - pos_0)
    if pos == 0 or pos == L:
        return 0
    else:
        return a * numpy.exp(b) * numpy.exp(c)


def V(x, x_start_step, V_0):
    """Returns potential for given x value."""
    if x >= x_start_step:
        return V_0
    else:
        return 0


def phi(pos, pos_0, delta_pos, E, sigma, x_start_step, V_0, prev_phi, delta_t):
    """ Calculates phi value """
    a = (1j*h_bar)/(2*m)
    b = psi(pos + delta_pos, pos_0, E, sigma) - 2 * psi(pos, pos_0, E, sigma) \
        + psi(pos - delta_pos, pos_0, E, sigma)
    c = delta_pos ** 2
    d = V(pos, x_start_step, V_0)
    e = psi(pos, pos_0, E, sigma)
    f = 1j * h_bar
    if pos == 0 or pos == L:
        return 0
    elif 0 < pos < L:
        return prev_phi + ((a * (b / c)) + ((d * e) / f)) * delta_t
    else:
        raise ValueError("Invalid position, must be between 0 and L.")


def prob_density(phi_value):
    """ Calculates probability density """
    return numpy.conj(phi_value) * phi_value


if __name__ == '__main__':
    sigma = 1E-8
    E = 0.2 * 1.602E-19
    V_0 = 0.16 * 1.602E-19
    x_start_step = 100E-9
    x_0 = 50E-9
    L = 200E-9

    delta_pos = 1.5E-10
    #delta_t = 2.25E-19
    #end_time = 5E-18
    time_steps = 2E6
    plot_step = 5000

    phi_values = []

    x_list = numpy.arange(0, L, delta_pos)

    #time_list = numpy.linspace(0, 40e-13, 100)
    #t = 0
    #while t < end_time:
    #    t_values.append(t)
    #    t += delta_t

    prob_densities = numpy.empty((len(time_list), len(x_list)))
    prev_phi = None
    x = 0
    for n in range(0, len(x_list)):
        for m in range(0, len(time_list)):
            if m == 0:
                current_phi = psi(x, x_0, E, sigma)
                prev_phi = current_phi
            else:
                x = x_list[n]
                current_phi = phi(
                    x, x_0, delta_pos, E, sigma,
                    x_start_step, V_0, prev_phi, delta_t)
                prev_phi = current_phi
            prob_densities.itemset((m, n), prob_density(current_phi))

    prob_density_values = pd.DataFrame(prob_densities, columns=x_list)
    V_values = [V(x, x_start_step, V_0) for x in x_list]

    for i in range(len(t_values)):
        plt.cla()
        plt.plot(prob_density_values.iloc[i])
        plt.plot(x_values, V_values)
        plt.show()
