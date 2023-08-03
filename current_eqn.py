import numpy as np
from parameters import *


def error_eqn_A5(phis_rad, phi_rad, alpha, Phi_ext_Phi0, epsilon_J):
    error_A5 = alpha * np.sin(phis_rad) + np.sin((phis_rad - 2 * np.pi * Phi_ext_Phi0) / 3) + epsilon_J * (
            M * phis_rad - phi_rad)
    return error_A5
