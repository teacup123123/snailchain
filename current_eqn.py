import numpy as np
import pylab as pl

from parameters import *
from scipy.optimize import fsolve


def residueI_eqnA5(phis_rad, phi_rad, alpha, Phi_ext_Phi0):
    """this reflects current conservation
    phis is phase across snail
    phi is phase across capacitor
    convention is current through M snails (each phis) and inductor (not free) share same direction
    capacitance is against them (phi)
    phiL = - M phis - phi
    """
    curr_alpha = alpha * np.sin(phis_rad)
    curr_unitaries = np.sin((phis_rad - 2 * np.pi * Phi_ext_Phi0) / 3)
    curr_snail = curr_unitaries + curr_alpha
    curr_inductor = Xi_J * (phi_rad - M * phis_rad)
    return curr_snail - curr_inductor  # should be zero


def fun_phis_of_phi(phi_rad, alpha, Phi_ext_Phi0):
    """given phi, phis will adapt itself to make residue equal to 0"""
    func = lambda x: residueI_eqnA5(x, phi_rad, alpha, Phi_ext_Phi0)
    phis = float(fsolve(func, [0.1]))
    return phis


if __name__ == '__main__':
    phi_rad_all = np.linspace(-100.*np.pi, 100.*np.pi, 101)
    phis_rad_all = [fun_phis_of_phi(phi_rad, alpha, Phi_ext_Phi0) for phi_rad in phi_rad_all]
    pl.plot(phi_rad_all, phis_rad_all)
    pl.show()
