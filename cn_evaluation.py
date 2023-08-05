import numpy as np
from scipy.optimize import fsolve
import math
from parameters import *
import U_Snail_M_arrays
import current_eqn
import pylab as pl
import differentiation
import U_SNAIL_as_a_func_phi_ext

#####################################################
""" This code is for evaluating the c1_tilde, c2_tilde, ...etc
coefficients corresponding to the potential energy term in Eq. A6 of Fratteni's article
Ref: https://arxiv.org/abs/1806.06093
"""


#####################################################

def c1_tilde_eq():
    [U_min_tilde, phi_min_rad_tilde] = U_SNAIL_as_a_func_phi_ext.phi_min_rad()
    phis_of_phi_min = current_eqn.fun_phis_of_phi(phi_min_rad_tilde, alpha=alpha, Phi_ext_Phi0=Phi_ext_Phi0)
    c1_tilde = Xi_J * (phi_min_rad_tilde - M * phis_of_phi_min)
    return c1_tilde


def c2_tilde():
    [U_min_tilde, phi_min_rad_tilde] = U_SNAIL_as_a_func_phi_ext.phi_min_rad()
    phis_of_phi_min = current_eqn.fun_phis_of_phi(phi_min_rad_tilde, alpha=alpha, Phi_ext_Phi0=Phi_ext_Phi0)
    phis_of_phi_func = lambda x: current_eqn.fun_phis_of_phi(x, alpha=alpha, Phi_ext_Phi0=Phi_ext_Phi0)
    First_der_phis_of_phi = differentiation.dern(phis_of_phi_func, 1, phi_min_rad_tilde)
    c2_tilde_val = Xi_J * (1 - M * First_der_phis_of_phi)
    return c2_tilde_val


def c3_tilde(Phi_ext_Phi0):
    [U_min_tilde, phi_min_rad_tilde] = U_SNAIL_as_a_func_phi_ext.phi_min_rad()
    phis_of_phi_min = current_eqn.fun_phis_of_phi(phi_min_rad_tilde, alpha=alpha, Phi_ext_Phi0=Phi_ext_Phi0)
    phis_of_phi_func = lambda x: current_eqn.fun_phis_of_phi(x, alpha=alpha, Phi_ext_Phi0=Phi_ext_Phi0)
    Sec_der_phis_of_phi = differentiation.dern(phis_of_phi_func, 2, phi_min_rad_tilde)
    c3_tilde_val = -M * Xi_J * Sec_der_phis_of_phi
    return c3_tilde_val


def c4_tilde(Phi_ext_Phi0):
    [U_min_tilde, phi_min_rad_tilde] = U_SNAIL_as_a_func_phi_ext.phi_min_rad()
    phis_of_phi_min = current_eqn.fun_phis_of_phi(phi_min_rad_tilde, alpha=alpha, Phi_ext_Phi0=Phi_ext_Phi0)
    phis_of_phi_func = lambda x: current_eqn.fun_phis_of_phi(x, alpha=alpha, Phi_ext_Phi0=Phi_ext_Phi0)
    Third_der_phis_of_phi = differentiation.dern(phis_of_phi_func, 3, phi_min_rad_tilde)
    c4_tilde_val = -M * Xi_J * Third_der_phis_of_phi
    return c4_tilde_val


if __name__ == '__main__':
    if 1:
        c1_tilde = c1_tilde_eq()
        print(c1_tilde)
        print(c2_tilde())
        # print(c3_tilde(Phi_ext_Phi0))
        #print(c4_tilde())

        Phi_ext_Phi0 = np.linspace(0,0.5,100)
        g3_tilde=[]
        g4_tilde = []
        for PeP0 in Phi_ext_Phi0:
            g3_GHz=phi_d1*phi_b1*phi_c2*c3_tilde(PeP0)*EJ_GHz
            g4_GHz=phi_d1*phi_b1*phi_c2**2*c4_tilde(PeP0)*EJ_GHz
            g3_tilde.append(abs(g3_GHz*1e3))
            g4_tilde.append(abs(g4_GHz*1e3))

        pl.plot(Phi_ext_Phi0,g3_tilde,'-o', label='M: {}'.format(M))
        pl.xlabel(r'$\Phi_{ext}/\Phi0$')
        pl.ylabel(r'$g_3$ (MHz)')
        pl.legend()
        pl.figure()
        pl.plot(Phi_ext_Phi0, g4_tilde, '-x', label='M: {}'.format(M))
        pl.xlabel(r'$\Phi_{ext}/\Phi0$')
        pl.ylabel(r'$g_4$ (MHz)')
        pl.legend()
        pl.show()
