# -*- coding: utf-8 -*-

__author__ = "Kjersti Rustad Kvisberg & Ida Lunde Naalsund"
__email__ = "kjkv@nmbu.no, idna@nmbu.no"

import numpy
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

def reflection_coef(x_values, phi_values):
    x_interval = int(len(x_values)/2)
    return abs(numpy.trapz(x_values[:x_interval],
                           numpy.conj(phi_values[:x_interval]) *
                                      phi_values[:x_interval]))

def transmission_coef(x_values, phi_values):
    x_interval = int(len(x_values)/2)
    return abs(numpy.trapz(x_values[x_interval:],
                           numpy.conj(phi_values[x_interval:]) *
                                      phi_values[x_interval:]))


if __name__ == '__main__':
    sigma = 1E-8
    E = 0.2 * 1.602E-19
    V_0 = 0.16 * 1.602E-19
    x_start_step = 100E-9
    x_0 = 50E-9
    L = 200E-9
    delta_pos = 1.5E-10
    delta_t = 2.25E-19
    time_steps = 2E6
    plot_step = 5000

    x_list = numpy.arange(0, L, delta_pos)
    psi_values = numpy.array([psi(pos, x_0, E, sigma) for pos in x_list])

    V_values = numpy.array([V(pos, x_start_step, V_0) for pos in x_list])
    a = delta_t / 1j * h_bar
    b = -h_bar ** 2 / 2 * m

    counter = 0
    img_num = 0

    for time in range(int(time_steps)):
        sec_deriv_psi =  (numpy.pad(psi_values[1:], (0, 1), 'constant',
                                constant_values=0) + numpy.pad
        (psi_values[0:-1],
                         (1, 0), 'constant', constant_values=0) - 2 *
                         psi_values) / delta_pos ** 2

        # Building the wave packet
        phi_values = psi_values + a * (b * sec_deriv_psi + V_values *
                                       psi_values)

        # Defining wave packet at boundaries
        phi_values[-1] = 0
        phi_values[0] = 0

        # Plotting the wave packet
        if counter % plot_step == 0:
            fig = plt.figure()
            plt.plot(x_list, (numpy.conj(phi_values)*phi_values))
            fig.savefig(f'img{str(img_num)}.png')
            plt.close(fig)
            plt.show()
            img_num += 1
            R = reflection_coef(x_list, phi_values)
            T = transmission_coef(x_list, phi_values)

        psi_values = phi_values
        counter += 1


