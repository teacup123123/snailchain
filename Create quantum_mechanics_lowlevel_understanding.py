import numpy as np
import pylab as pl


def psi_zero(x):
    return 0.


def psi(x):
    return np.exp(-x ** 2 / 2) * np.exp(1j * x)


def psi0(x):
    return np.exp(-x ** 2 / 2)


def psi1(x):
    return np.exp(-x ** 2 / 2) * x


def psi2(x):
    return np.exp(-x ** 2 / 2) * np.exp(1j * x) * np.cos(x ** 2)


def psi3(x):
    return np.exp(-x ** 2 / 2) * np.exp(1j * x) * np.cos(x ** 3)


def op_trans_by_2pi(psi):
    "S_2pi = f->(x->f(x+2pi))"
    # def translated_psi(x):
    #     return psi(x + 2 * np.pi)
    return lambda x: psi(x + 2 * np.pi)


def op_flipx(psi):
    "f->(x->f(-x))"
    return lambda x: psi(-x)


epsilon = 0.001


def Hho(psi):
    "f->(x->(x^2-dx^2)f(x))"

    def Hf(x):
        first_part = x * x * psi(x)
        second_part = (psi(x + epsilon) - 2 * psi(x) + psi(x - epsilon)) / epsilon ** 2
        return first_part - second_part

    return Hf


def Uho(psi):
    "f->(x->(x^2)f(x))"

    def Hf(x):
        first_part = x * x * psi(x)
        return first_part

    return Hf


def proba_amplitude(psi):
    def pa(x):
        return np.abs(psi(x)) ** 2

    return pa


def psi_cosphi(x):
    return np.cos(x)


def U_minusEjcosphi(psi):
    # equivalent to modulate_by(psi, lambda x:-Ej*np.cos(x))
    Ej = 1
    return lambda x: -Ej * np.cos(x) * psi(x)


def modulate_by(psi, f):
    return lambda x: f(x) * psi(x)


def Ucos(psi):
    return modulate_by(psi, np.cos)

def commutator(op1, op2):
    def returned(psi):
        return lambda x:op1(op2(psi))(x)-op2(op1(psi))(x)
    return returned

def commutator_Hho_opflipx(psi):
    #equivalent to commutator(Hho, op_flipx)
    return lambda x: Hho(op_flipx(psi))(x) - op_flipx(Hho(psi))(x)


def main_continuous():
    x = np.linspace(-10., 10., 201) # rough approx of entire real space
    # tests =[psi, op_flipx(psi), op_trans_by_2pi(psi), commutator_Hho_opflipx(psi)]
    tests = [psi1, U_minusEjcosphi(psi1)] # to be visualized
    for wavefun_2b_plotted in tests:
        pl.figure()
        pl.plot(x, wavefun_2b_plotted(x).real)
        pl.plot(x, wavefun_2b_plotted(x).imag)
    pl.show()

    tests = [psi, psi2, psi3, psi0]  # .... rough approximation of all possible wavefunctions
    coords = [0., 4., np.pi]  # .... # rough approximation of the entire space
    truth = all((all(commutator_Hho_opflipx(f)(x) == psi_zero(x) for x in coords)) for f in tests) # this says [Hho, opflipx] = 0 (0 being the operator 0)
    print(truth)
    print()

    # pl.plot(x, proba_amplitude(psi)(x))


def operator_p(psi):
    return lambda x: -1j * (psi(x + epsilon) - psi(x - epsilon)) / (2 * epsilon)


def psik(x, k):
    # wave function with momentum k
    # assert -np.pi <= x <= np.pi
    return np.exp(1j * k * x)


def superconductors():
    coords = x = np.arange(-np.pi, np.pi, np.pi / 25)

    tests = [
        lambda x: psik(x, 3),
        # modulate_by((lambda x: psik(x, 3)), np.cos),
        Ucos(lambda x: psik(x, 3)),
        lambda x: psik(x, 2) / 2 + psik(x, 4) / 2,
    ]
    for wavefun_2b_plotted in tests:
        pl.figure()
        pl.plot(x, wavefun_2b_plotted(x).real)
        pl.plot(x, wavefun_2b_plotted(x).imag)
    pl.show()


if __name__ == '__main__':
    superconductors()
    # assert False
    print('program doing fine')
    # main_continuous()
