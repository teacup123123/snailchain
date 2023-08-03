import numpy as np
from parameters import *


def residueI_eqnA5(phis_rad, phi_rad, alpha, Phi_ext_Phi0, epsilon_J):
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
    curr_inductor = phi_rad - M * phis_rad
    return curr_snail - curr_inductor  # should be zero
