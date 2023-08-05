import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import minimize
import math
from parameters import *
import U_Snail_M_arrays
import current_eqn
import pylab as pl

###############################################################################
""" This is to plot U_SNAIL as a function of Phi_ext to get
    an intuitive feeling of the function minima and maxima points. 
    Also, to find out the phi_min_rad minimum corresponding to the potential. 
    """


###############################################################################
def U_of_phi(phi_rad, Phi_ext_Phi0, alpha, ns, M, EL, EJ):
    """
    This is to plot U_SNAIL as a function of Phi_ext to get
    an intuitive feeling of the function minima and maxima points
    """
    func = lambda x: current_eqn.residueI_eqnA5(x, phi_rad, alpha, Phi_ext_Phi0)
    phis = fsolve(func, [0.0])[0]
    U_snail = U_Snail_M_arrays.U_SNAIL_by_EJ(phis, phi_rad, Phi_ext_Phi0, alpha, ns, M, EL, EJ)
    return U_snail


def phi_min_rad():
    """ To evaluate the phi_min_rad correponding to the SNAIL potential"""
    phi_min_rad_result = minimize(lambda x: U_of_phi(x, Phi_ext_Phi0, alpha, ns, M, EL, EJ),
                                  x0=3.5)  # x0 is the initial guess for the minimum
    return [phi_min_rad_result.fun, phi_min_rad_result.x[0]]


if __name__ == '__main__':
    phi_rad_all = np.linspace(-15. * np.pi, 15. * np.pi, 101)
    if 0:
        Phi_ext_Phi0 = 2 * np.pi * np.linspace(-2 * np.pi, 2 * np.pi, 101)
        data = np.zeros((len(phi_rad_all), len(Phi_ext_Phi0)))
        for n_PeP0, PeP0 in enumerate(Phi_ext_Phi0):
            for n_phi, phi in enumerate(phi_rad_all):
                U_s = U_of_Phi_ext(phi, PeP0, alpha, ns, M, EL, EJ)

        # Create the pcolor plot
        pl.pcolormesh(Phi_ext_Phi0, phi_rad_all, data, cmap='viridis')
        pl.colorbar(label='Value')
        pl.xlabel('X')
        pl.ylabel('Y')
        pl.title('Potential energy')
        pl.show()

    if 1:
        U_s = [U_of_phi(phi, Phi_ext_Phi0, alpha, ns, M, EL, EJ) for phi in zip(phi_rad_all)]
        # Find the minimum using minimize
        [U_min_val, phi_rad_min_val] = phi_min_rad()

        # Display the result
        print("Minimum value of the function:", U_min_val)
        print("phi_rad-coordinate of the minimum:", phi_rad_min_val)

        pl.plot(phi_rad_all, U_s)
        # Mark the minimum on the plot
        pl.plot(phi_rad_min_val, U_min_val, 'ro', label='Minimum')
        pl.show()
