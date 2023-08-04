import numpy as np
import pylab as pl

from parameters import *
from scipy.optimize import fsolve


def curr_snail(phis_rad, alpha=alpha, Phi_ext_Phi0=Phi_ext_Phi0):
    curr_alpha = alpha * np.sin(phis_rad)
    curr_unitaries = np.sin((phis_rad - 2 * np.pi * Phi_ext_Phi0) / 3)
    curr_snail = curr_unitaries + curr_alpha
    return curr_snail


# def der_curr_snail(phis_rad, alpha=alpha, Phi_ext_Phi0=Phi_ext_Phi0):
#     return alpha * np.cos(phis_rad) + np.cos((phis_rad - 2 * np.pi * Phi_ext_Phi0) / 3) / 3


def curr_inductor(phis_rad, phi_rad):
    # - Xi_J * M * phis_rad + Xi_J phi_rad
    return Xi_J * (phi_rad - M * phis_rad)


# def der_curr_inductor():
#     return -Xi_J * M


def residueI_eqnA5(phis_rad, phi_rad, alpha=alpha, Phi_ext_Phi0=Phi_ext_Phi0):
    """this reflects current conservation
    phis is phase across snail
    phi is phase across capacitor
    convention is current through M snails (each phis) and inductor (not free) share same direction
    capacitance is against them (phi)
    phiL = - M phis - phi
    """
    print('eval')
    return curr_snail(phis_rad, alpha, Phi_ext_Phi0) - curr_inductor(phis_rad, phi_rad)  # should be zero


# def der_residueI_eqnA5(phis_rad, alpha=alpha, Phi_ext_Phi0=Phi_ext_Phi0):
#     return 2 * (der_curr_snail(phis_rad, alpha, Phi_ext_Phi0) - der_curr_inductor())  # should be zero


def fun_phis_of_phi(phi_rad, alpha=alpha, Phi_ext_Phi0=Phi_ext_Phi0):
    """given phi, phis will adapt itself to make residue equal to 0"""
    func = lambda x: residueI_eqnA5(x, phi_rad, alpha, Phi_ext_Phi0)
    # phis = float(fsolve(func, [0.0], fprime=lambda x: der_residueI_eqnA5(x, alpha, Phi_ext_Phi0)))
    phis = float(fsolve(func, [0.0]))
    return phis


if __name__ == '__main__':
    phis_all = np.linspace(-2. * np.pi, 2. * np.pi, 201)
    pl.plot(phis_all, curr_snail(phis_all))
    pl.plot(phis_all, curr_inductor(phis_all, -0.01), label=-0.01)
    pl.plot(phis_all, curr_inductor(phis_all, 0.0), label=0.)
    pl.plot(phis_all, curr_inductor(phis_all, 0.01), label=0.01)
    pl.plot(phis_all, curr_inductor(phis_all, 0.02), label=0.02)
    pl.plot(phis_all, curr_inductor(phis_all, 2 * np.pi), label=2 * np.pi)
    pl.title([fun_phis_of_phi(x) for x in [-.01, 0, .01, .02, 2 * np.pi]])
    pl.xlabel('phis')
    pl.ylim([-1, 1])
    pl.legend()
    pl.show()

    if 0:
        phis_rad_1 = np.linspace(-1 * np.pi, 1 * np.pi, 101)
        phi_rad_1 = 2 * np.pi
        curr_alpha_1 = alpha * np.sin(phis_rad_1)
        curr_unitaries = np.sin((phis_rad_1 - 2 * np.pi * Phi_ext_Phi0) / 3)
        curr_snail_1 = curr_alpha_1 + curr_unitaries
        curr_inductor = Xi_J * (phi_rad_1 - M * phis_rad_1)

        pl.plot(phis_rad_1, curr_alpha_1, label='I_snail', color='blue', linestyle='-', marker='o')
        pl.plot(phis_rad_1, curr_inductor, label='I_inductor', color='green', linestyle='--', marker='x')
        pl.title('Multiple Graphs in One Plot')
        pl.xlabel('X-axis')
        pl.ylabel('Y-axis')
        pl.legend()
        pl.show()

    if 0:
        phis_rad_all = [fun_phis_of_phi(phi_rad, alpha, Phi_ext_Phi0) for phi_rad in phi_rad_all]
        pl.plot(phi_rad_all, phis_rad_all)
        pl.show()
