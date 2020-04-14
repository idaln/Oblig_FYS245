# -*- coding: utf-8 -*-

__author__ = "Ida Lunde Naalsund & Kjersti Rustad Kvisberg"
__email__ = "idna@nmbu.no, kjkv@nmbu.no"
import math
import numpy as np
import matplotlib.pyplot as plt


# Oppgave2
sigma = 1 * 10 ** -8
E = 0.2 * 1.60 * 10 ** -19
V_0 = 0.16 * 1.60 * 10 ** -19
m = 9.11 * 10 ** -31
x_0 = 50 * 10 ** -9
L = 200 * 10 ** -9
delta_x = 1.5 * 10 ** -10
delta_t = 2.25 * 10 ** -19
step_start = 100 * 10 ** -9
h_bar = 6.62607 * 10 ** -34 / (2 * math.pi)

timesteps = 2 * 10 ** 6
plot_step = 5000


def k(x):
    return (1 / h_bar) * math.sqrt(2 * m * (E - V(x)))


def psi(x):
    a = np.exp((-(x - x_0) ** 2) / (2 * sigma ** 2))
    b = np.exp(1j * k(x) * (x - x_0))
    return (1 / ((np.pi ** 0.25) * np.sqrt(sigma))) * a * b


def V(x):
    if 0 <= x < step_start:
        return 0
    else:
        return V_0

# Lager video
#import cv2
#import os
#
# image_folder = 'C:/Users/idaln/OneDrive/Dokumenter/FYS245'
#
# video_name = 'Kvantepropagering.avi'
#
# images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
# frame = cv2.imread(os.path.join(image_folder, images[0]))
# height, width, layers = frame.shape
#
# video = cv2.VideoWriter(video_name, 0, 1, (width, height))
#
# for image in images:
#     video.write(cv2.imread(os.path.join(image_folder, image)))
#
# video.release()
# cv2.destroyAllWindows()

if __name__ ==  "__main__":

    x_v = np.arange(0, L, delta_x)
    psi_v = np.array([psi(x) for x in x_v])
    psi_v[len(x_v) - 1] = 0
    psi_v[0] = 0
    V_values = np.array([V(x) for x in x_v])
    a = (delta_t / (1j * h_bar))
    b = (-(h_bar ** 2 / (2 * m)))
    counter = 0
    pict_count = 0
    for time in range(timesteps):
        d2psi = (np.pad(psi_v[1:len(x_v)], (0, 1), 'constant',
                        constant_values=0) + np.pad(psi_v[0:len(x_v) - 1],
                                                    (1, 0),
                                                    'constant',
                                                    constant_values=0) - 2 * psi_v) / delta_x ** 2
        V_psi = V_values * psi_v
        Psi_v = psi_v + a * (b * d2psi + V_psi)
        Psi_v[len(x_v) - 1] = 0
        Psi_v[0] = 0
        if counter % plot_step == 0:
            fig = plt.figure()
            plt.plot(x_v, Psi_v * np.conj(Psi_v))
            fig.savefig(f"img{str(pict_count)}.png")
            plt.close(fig)
            plt.show()
            pict_count += 1
            x1 = x_v[:int(len(x_v) / 2)]
            x2 = x_v[int(len(x_v) / 2):]
            R_coeff = abs(np.trapz(x1, Psi_v[:int(len(x_v) / 2)] * np.conj(
                Psi_v[:int(len(x_v) / 2)])))
            T_coeff = abs(np.trapz(x2, Psi_v[int(len(x_v) / 2):] * np.conj(
                Psi_v[int(len(x_v) / 2):])))
        psi_v = Psi_v
        counter += 1

#def phi(pos, pos_0, delta_pos, E, sigma, x_start_step, V_0, prev_phi, delta_t):
#    """ Calculates phi value """
#    a = (1j*h_bar)/(2*m)
#    b = psi(pos + delta_pos, pos_0, E, sigma) - 2 * psi(pos, pos_0, E, sigma) \
#        + psi(pos - delta_pos, pos_0, E, sigma)
#    c = delta_pos ** 2
#    d = V(pos, x_start_step, V_0)
#    e = psi(pos, pos_0, E, sigma)
#    f = 1j * h_bar
#    if pos == 0 or pos == L:
#        return 0
#    elif 0 < pos < L:
#        return prev_phi + ((a * (b / c)) + ((d * e) / f)) * delta_t
#    else:
#        raise ValueError("Invalid position, must be between 0 and L.")


#def prob_density(phi_value):
#    """ Calculates probability density """
#    return numpy.conj(phi_value) * phi_value

#    prob_densities = numpy.empty((len(time_list), len(x_list)))
#    prev_phi = None
#    x = 0
#    for n in range(0, len(x_list)):
#        for m in range(0, len(time_list)):
#            if m == 0:
#                current_phi = psi(x, x_0, E, sigma)
#                prev_phi = current_phi
#            else:
#                x = x_list[n]
#                current_phi = phi(
#                    x, x_0, delta_pos, E, sigma,
#                    x_start_step, V_0, prev_phi, delta_t)
#                prev_phi = current_phi
#            prob_densities.itemset((m, n), prob_density(current_phi))

#    prob_density_values = pd.DataFrame(prob_densities, columns=x_list)
#    V_values = [V(x, x_start_step, V_0) for x in x_list]

