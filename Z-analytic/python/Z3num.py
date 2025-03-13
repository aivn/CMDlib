"""
Численный расчет трехчастичной статистической суммы.
"""
from Z3backend import *
import Z2num
import numpy as np
from scipy.integrate import quad
import scipy.optimize as opt

_z3 = Z3_symmetrical()


def calc_from_coeffs(h_i, h_j, l_ij, l_ik):
    """Функциональная обертка над Z3_symmetrical.calc_from_coeffs. Возвращает объект класса Z3_symmetrical."""
    _z3.h_i, _z3.h_j, _z3.l_ij, _z3.l_ik = h_i, h_j, l_ij, l_ik
    _z3.calc_from_coeffs()
    return _z3


def coeffs_symmetrical_approx(m, eta_ij, eta_ik):
    """Аппроксимация коэффициентов для симметричного случая."""
    if m == 1 or eta_ij == 1 or eta_ik == 1:
        return np.nan, np.nan, np.nan, np.nan

    zeta_ij = (eta_ij - m * m) / (1 - m * m)
    zeta_ik = (eta_ik - m * m) / (1 - m * m)

    zeta_1 = zeta_ij
    zeta_2 = (zeta_ik - zeta_ij * zeta_ij) / (1 - zeta_ij * zeta_ij)

    rho_h_i = 0.92356815 - 0.0079533 * m  - 0.00436018 * zeta_1 - 0.46472716 * zeta_2
    rho_h_j = 0.88135978 - 0.00869924 * m - 0.80522817 * zeta_1 + 0.92599076 * zeta_2
    rho_l_ij = 9.20503066e-01 - 6.08982188e-04 * m - 1.14971741e-02 * zeta_1 - 4.62352724e-01 * zeta_2
    rho_l_ik = 0.69791727 + 0.01027317 * m - 0.37658349 * zeta_1 + 0.57813 * zeta_2

    Z2_h_ij, Z2_l_ij = Z2num.coeffs_symmetrical_approx(m, eta_ij)
    Z2_h_ik, Z2_l_ik = Z2num.coeffs_symmetrical_approx(m, eta_ik)

    return rho_h_i * Z2_h_ij, rho_h_j * Z2_h_ik, rho_l_ij * Z2_l_ij, rho_l_ik * Z2_l_ik


def find_coeffs_symmetrical(m, eta_ij, eta_ik, eps=1e-8, max_steps=64):
    """Функциональная обертка над Z3_symmetrical.calc_from_moments. Дополнительно использует аппроксимации коэффициентов. Возвращает объект класса Z3_symmetrical."""
    _z3.mi, _z3.mj, _z3.eta_ij, _z3.eta_ik = m, m, eta_ij, eta_ik
    _z3.calc_from_moments(*coeffs_symmetrical_approx(m, eta_ij, eta_ik), eps=eps, max_steps=max_steps)
    return _z3


def find_coeffs_symmetrical_by_scipy(m, eta_ij, eta_ik):
    """См. find_coeffs_symmetrical. Использует scipy для решения обратной задачи."""
    def func(args):
        _z3.h_i, _z3.h_j, _z3.l_ij, _z3.l_ik = args
        _z3.calc_from_coeffs()
        return [_z3.mi-m, _z3.mj-m, _z3.eta_ij-eta_ij, _z3.eta_ik-eta_ik]

    result = opt.fsolve(func, coeffs_symmetrical_approx(m, eta_ij, eta_ik), full_output=True)
    if result[2] == 1:
        coeffs = result[0]
    else:
        coeffs = np.nan, np.nan, np.nan, np.nan

    _z3.h_i, _z3.h_j, _z3.l_ij, _z3.l_ik = coeffs
    _z3.calc_from_coeffs()
    return _z3

# =======================================================================================

def calc_Z3_integrate(hi, hj, hk, lij, ljk, lik):
    def abs2_i(x): return hi**2+lij**2+2*hi*lij*x
    def abs2_k(x): return hk**2+ljk**2+2*hk*ljk*x

    def dot_i_k(x): return hi*hk+lij*ljk+x*(hi*ljk+hk*lij)

    def sqrt_i_z(x): return dot_i_k(x)/np.sqrt(abs2_k(x))

    def sqrt_i_y(x):
        if hi == hk and lij == ljk:
            return 0
        return np.sqrt(abs2_i(x) - dot_i_k(x)**2/abs2_k(x))

    def func(x): return np.exp(hj*x)*Z2num.calc_Z2_integrate(sqrt_i_y(x), sqrt_i_z(x), np.sqrt(abs2_k(x)), lik)
    return 2*np.pi*quad(func, -1, 1)[0]
