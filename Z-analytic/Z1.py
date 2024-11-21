import numpy as np
import scipy.optimize as opt


def mu_p_approx(m):
    m2 = m * m
    mu_p_arg = m * ((0.775383) + ((0.109185) + ((0.289114) + ((-0.871214) +
                    ((1.85968) + ((-1.71306) + (0.550912) * m2) * m2) * m2) * m2) * m2) * m2)
    mu_p = (1 - mu_p_arg * mu_p_arg) / 3
    return mu_p


def calc_Z1(p):
    return 4*np.pi*np.sinh(p)/p


def calc_m(p):
    return 1/np.tanh(p)-1/p


def find_coeffs(m):
    def func(p):
        return calc_m(p)-m

    mu_p = mu_p_approx(m)
    p = m / mu_p

    result = opt.root(func, p)
    if result["success"]:
        result = result["x"]
    else:
        result = np.asarray([np.nan])
    return result


find_coeffs = np.vectorize(find_coeffs, signature='()->()')
