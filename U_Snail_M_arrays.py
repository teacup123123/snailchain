import numpy as np
from parameters import *


def U_SNAIL_by_EJ(phis_rad, phi_rad, Phi_ext_Phi0, alpha, ns, M, EL, EJ):
    """ This Eq. is the potential energy term corresponding to
    SNAIL arrays + inductor L
    """
    U_single_snail_by_EJ = - alpha * np.cos(phis_rad) - 3 * np.cos((2 * np.pi * Phi_ext_Phi0 - phis_rad) / ns)
    U_inductor_by_EJ = (0.5 * Xi_J) * (phi_rad - M * phis_rad) ** 2
    U_SPA_by_EJ_Eq = M * U_single_snail_by_EJ + U_inductor_by_EJ
    return U_SPA_by_EJ_Eq
