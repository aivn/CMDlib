import numpy as np
import numba as nb


def single_and_vectorize(func):
    return nb.njit(func, cache=True), nb.vectorize([nb.float64(*[nb.float64]*func.__code__.co_argcount)], target="parallel", cache=True)(func)


# ======================================== MOMENTS ========================================
def calc_Z1_norm(p):
    if p == 0:
        return 0
    return 2*np.pi*(1-np.exp(-2*p))/p


def calc_m(p):
    if p == 0:
        return 0
    return 1/np.tanh(p)-1/p


calc_Z1_norm_single, calc_Z1_norm = single_and_vectorize(calc_Z1_norm)
calc_m_single, calc_m = single_and_vectorize(calc_m)


# ======================================== FINDING ========================================
@nb.njit(cache=True)
def mu_p_approx(m):
    m2 = m * m
    mu_p_arg = m * ((0.775383) + ((0.109185) + ((0.289114) + ((-0.871214) + ((1.85968) + ((-1.71306) + (0.550912) * m2) * m2) * m2) * m2) * m2) * m2)
    mu_p = (1 - mu_p_arg * mu_p_arg) / 3
    return mu_p


@nb.njit(cache=True)
def coeffs_approx(m):
    return m / mu_p_approx(m)


@nb.njit(cache=True)
def find_coeffs_single(m, eps=1e-8):
    if m == 0:
        return 0

    if m == 1:
        return np.nan

    p = coeffs_approx(m)
    delta_p = 1

    while abs(delta_p) > eps:
        _exp_p = np.exp(-2*p)
        coth = (1+_exp_p)/(1-_exp_p)

        _p = 1/p
        _p2 = _p*_p
        coth2 = coth*coth

        expr_m = coth - _p
        expr_m_diff_p = 1 - coth2 + _p2

        delta_p = (expr_m - m) / expr_m_diff_p
        p -= delta_p

    return p


@nb.njit(parallel=True, cache=True)
def find_coeffs_vectorize(m, eps=1e-8):
    p = np.empty_like(m)

    for i in nb.prange(m.size):
        p.flat[i] = find_coeffs_single(m.flat[i], eps)

    return p


def find_coeffs(m, eps=1e-8):
    if isinstance(m, np.ndarray):
        return find_coeffs_vectorize(m, eps)
    else:
        return find_coeffs_single(m, eps)
