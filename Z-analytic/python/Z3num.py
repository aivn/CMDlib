"""
Численный расчет трехчастичной статистической суммы.
"""
from Z3backend import *
import Z2num
import numpy as np
from scipy.integrate import quad


def calc_Z3_integrate(h_i, h_j, l_ij, l_ik):
    _sqrt = lambda x: np.sqrt(h_i**2+l_ij**2+2*h_i*l_ij*x)
    func = lambda x: np.exp(h_j*x)*Z2num.calc_Z2_integrate(0, _sqrt(x), _sqrt(x), l_ik)
    return 2*np.pi*quad(func, -1, 1)[0]