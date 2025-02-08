"""
Численный расчет двухчастичной статистической суммы.
"""
from Z2backend import cF2_norm_symmetrical, E_norm_symmetrical, Z2_symmetrical
import scipy.optimize as opt
import numpy as np

_z2 = Z2_symmetrical()


def coeffs_symmetrical_approx(m, eta):
    """Аппроксимация коэффициентов для симметричного случая."""
    m2 = m * m
    zeta = (eta - m2) / (1 - m2)
    zeta2 = zeta * zeta

    mu_p_arg = m * ((0.775383) + ((0.109185) + ((0.289114) + ((-0.871214) +
                    ((1.85968) + ((-1.71306) + (0.550912) * m2) * m2) * m2) * m2) * m2) * m2)
    mu_p = (1 - mu_p_arg * mu_p_arg) / 3

    upsilonM = (1 - eta) / (1 - m2) * m * mu_p + 0.0467246 * m * (1 - m) * zeta * (1 - zeta) * (1 - 7.96819 * m - 0.939775 * zeta + 14.0242 * m * zeta + 20.9323 * m2 + 7.65054 * zeta2 - 14.5937 *
                                                                                                m2 * zeta - 5.35767 * m * zeta2 - 13.542 * m * m2 - 3.72091 * zeta * zeta2) + (0.00349416 + (0.0191611 + (-0.0367299 + (-0.0684992 + (0.149972 - 0.067404 * m2) * m2) * m2) * m2) * m2) * m
    upsilon = upsilonM / m
    rho = 1 - .5 * zeta - .514 * (1 - zeta) * zeta * (1 - (.232 + (1.681 - 1.466 * m) * m)
                                                      * m + (-.123 * m - 1.041 + (0.858 - 0.356 * zeta) * zeta) * zeta)

    return rho * m / mu_p, (1 - rho) / upsilon


def calc_from_coeffs_symmetrical(h, l):
    """Функциональная обертка над Z2_symmetrical.calc_from_coeffs. Возвращает объект класса Z2_symmetrical."""
    _z2.h, _z2.l = h, l
    _z2.calc_from_coeffs()
    return _z2


def find_coeffs_symmetrical(m, eta):
    """Функциональная обертка над Z2_symmetrical.calc_from_moments. Дополнительно использует аппроксимации коэффициентов. Возвращает объект класса Z2_symmetrical."""
    _z2.m, _z2.eta = m, eta
    _z2.calc_from_moments(*coeffs_symmetrical_approx(m, eta))
    return _z2


def find_coeffs_symmetrical_by_scipy(m, eta):
    """См. find_coeffs_symmetrical. Использует scipy для решения обратной задачи."""
    def func(args):
        _z2.h, _z2.l = args
        _z2.calc_from_coeffs()
        return [_z2.m-m, _z2.eta-eta]

    result = opt.fsolve(func, coeffs_symmetrical_approx(m, eta), full_output=True)
    if result[2] == 1:
        coeffs = result[0]
    else:
        coeffs = np.nan, np.nan

    _z2.h, _z2.l = coeffs
    _z2.calc_from_coeffs()
    return _z2
