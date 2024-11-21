from scipy.special import dawsn
import numpy as np
from functools import lru_cache
import scipy.optimize as opt
import Z1


SIGMAS = [-1, 1]


def F2(x):
    return dawsn(x)


@lru_cache(maxsize=1)
def _Z2_E_terms_norm(h, l):
    def func(s1, s2): return np.exp((s1+s2)*h+s1*s2*l-2*h-l)
    return np.asarray([[func(s1, s2) for s2 in SIGMAS] for s1 in SIGMAS])


@lru_cache(maxsize=1)
def _Z2_F_terms(h, l):
    def func(s1, s2): return F2((h+(s1+s2)*l)/np.sqrt(2*l))
    return np.asarray([[func(s1, s2) for s2 in SIGMAS] for s1 in SIGMAS])


def _Z2_sum_terms(terms, sigmas_func):
    result = 0
    for s1 in range(2):
        for s2 in range(2):
            result += sigmas_func(SIGMAS[s1], SIGMAS[s2])*terms[s1, s2]
    return result


@lru_cache(maxsize=1)
def f_cF2Function_norm(h, l):
    _terms = _Z2_E_terms_norm(h, l)*_Z2_F_terms(h, l)
    def _sigmas_func_12(s1, s2): return s1*s2
    return _Z2_sum_terms(_terms, _sigmas_func_12)


@lru_cache(maxsize=1)
def EFunction_norm(h, l, arg_s1, arg_s2):
    def _sigmas_func(s1, s2):
        result = 1
        if arg_s1:
            result *= s1
        if arg_s2:
            result *= s2
        return result
    return _Z2_sum_terms(_Z2_E_terms_norm(h, l), _sigmas_func)


def calc_Z2_norm(h, l):
    return np.sqrt(2)*(2*np.pi)**2/(h*np.sqrt(l))*f_cF2Function_norm(h, l)


def calc_m(h, l):
    return -1/(2*h) - h/(2*l) + np.sqrt(2)*EFunction_norm(h, l, 1, 1)/(4*np.sqrt(l)*f_cF2Function_norm(h, l))


def calc_eta(h, l):
    return -1 - 1/(2*l) + h**2/(2*l**2) + np.sqrt(2)*EFunction_norm(h, l, 1, 0)/(2*np.sqrt(l)*f_cF2Function_norm(h, l)) - np.sqrt(2)*h*EFunction_norm(h, l, 1, 1)/(4*l**(3/2)*f_cF2Function_norm(h, l))


def find_coeffs(m, eta):
    if m == 0 or eta == 0:
        return np.asarray([np.nan, np.nan])

    if m == 1 or eta == 1:
        return np.asarray([np.nan, np.nan])

    def func(coeffs):
        return [
            calc_m(*coeffs)-m,
            calc_eta(*coeffs)-eta
        ]

    m2 = m * m
    zeta = (eta - m2) / (1 - m2)
    zeta2 = zeta * zeta

    mu_p = Z1.mu_p_approx(m)

    upsilonM = (1 - eta) / (1 - m2) * m * mu_p + 0.0467246 * m * (1 - m) * zeta * (1 - zeta) * (1 - 7.96819 * m - 0.939775 * zeta + 14.0242 * m * zeta + 20.9323 * m2 + 7.65054 * zeta2 - 14.5937 * m2 * zeta - 5.35767 * m * zeta2 - 13.542 * m * m2 - 3.72091 * zeta * zeta2) + (0.00349416 + (0.0191611 + (-0.0367299 + (-0.0684992 + (0.149972 - 0.067404 * m2) * m2) * m2) * m2) * m2) * m
    upsilon = upsilonM / m
    rho = 1 - .5 * zeta - .514 * (1 - zeta) * zeta * (1 - (.232 + (1.681 - 1.466 * m) * m) * m + (-.123 * m - 1.041 + (0.858 - 0.356 * zeta) * zeta) * zeta)

    h = rho * m / mu_p
    l = (1 - rho) / upsilon

    result = opt.root(func, [h, l])
    if result["success"]:
        result = result["x"]
    else:
        result = np.asarray([np.nan, np.nan])
    return result


calc_Z2_norm = np.vectorize(calc_Z2_norm)
calc_m = np.vectorize(calc_m)
calc_eta = np.vectorize(calc_eta)
find_coeffs = np.vectorize(find_coeffs, signature='(),()->(n)')
