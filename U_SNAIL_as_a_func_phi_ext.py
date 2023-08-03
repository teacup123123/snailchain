import numpy as np
from scipy.optimize import fsolve
import math
from parameters import *
import U_Snail_M_arrays
import current_eqn
import pylab as pl


def U_of_Phi_ext(phi_rad, Phi_ext_Phi0, alpha, ns, M, EL, EJ):
    """
    This is to plot U_SNAIL as a function of Phi_ext to get
    an intuitive feeling of the function minima and maxima points
    """
    func = lambda x: current_eqn.residueI_eqnA5(x, phi_rad, alpha, Phi_ext_Phi0)
    phis = fsolve(func, [0.1])[0]
    U_snail = U_Snail_M_arrays.U_SNAIL_by_EJ(phis, phi_rad, Phi_ext_Phi0, alpha, ns, M, EL, EJ)
    return U_snail


if __name__ == '__main__':
    Phi_ext_Phi0 = 2 * np.pi * np.linspace(-2 * np.pi, 2 * np.pi, 101)
    phi_rad_all = np.linspace(-1000. * np.pi, 1000. * np.pi, 101)
    data = np.zeros((len(phi_rad_all), len(Phi_ext_Phi0)))
    for PeP0, n_PeP0 in enumerate(Phi_ext_Phi0):
        for phi, n_phi in enumerate(phi_rad_all):
            U_s = U_of_Phi_ext(phi, PeP0, alpha, ns, M, EL, EJ)
    # pl.plot(phi_rad_all, U_s)
    # pl.show()

    # Create the pcolor plot
    pl.pcolormesh(Phi_ext_Phi0, phi_rad_all, data, cmap='viridis')
    pl.colorbar(label='Value')
    pl.xlabel('X')
    pl.ylabel('Y')
    pl.title('Pcolor Plot with Two Variable For Loops')
    pl.show()
