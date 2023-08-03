############ system params ###########################
# SNAIL
from constants import *

M, alpha, ns = 31, 0.102, 3
LJ = 16.86e-9
EJ_Hz = phi0 ** 2 / LJ
EJ = h * EJ_Hz
# Cavity
w_cavity = 2 * np.pi * 13e9
C_cavity = 9e-13  # obtained from an approximate expression for C=epsilon0*pi*R0^2/d for reentrant cavity
L_cavity = 1 / (w_cavity ** 2 * C_cavity)
EL = 4 * np.pi ** 2 * phi0 ** 2 / L_cavity  # is it right?
# Common
Xi_J = LJ / L_cavity
Phi_ext_Phi0 = 2 * np.pi * 0.25  # np.linspace(0, 0.5, 100)

# HFSS and pyEPR sim results
fd_Hz, fb_Hz, fc_Hz = 5.92639e9, 6.29573e9, 12.0198e9
participation_matrix = np.array([[0.00088, 0.45, 0.52], [8e-7, 0.53, 0.45], [2.4e-10, 0.00055, 0.00056]])
freq_Hz = np.array([[fb_Hz, 0, 0], [0, fd_Hz, 0], [0, 0, fc_Hz]])
w_radHz = 2 * np.pi * freq_Hz
const = np.sqrt(hbar / (2 * EJ_Hz * h))

# ZPFs (these values may change)
phi_d1, phi_b1, phi_c2 = 0.3704, 0.4143, 0.0186

###########################################################
