from scipy.special import dawsn  # requires numba-scipy
import numpy as np
import numba as nb
import Z1

# ======================================== BASE ========================================
SIGMAS = np.asarray([-1, 1])


@nb.njit
def F2(x):
    return dawsn(x)


@nb.njit
def _Z2_E_terms_norm(h, l):
    def func(s1, s2): return np.exp((s1+s2)*h+s1*s2*l-2*h-l)
    return np.asarray([[func(s1, s2) for s2 in SIGMAS] for s1 in SIGMAS])


@nb.njit
def _Z2_F_terms(h, l):
    def func(s1, s2): return F2((h+(s1+s2)*l)/np.sqrt(2*l))
    return np.asarray([[func(s1, s2) for s2 in SIGMAS] for s1 in SIGMAS])


@nb.njit
def _Z2_sum_terms(terms, sigmas_func, *sigmas_args):
    result = 0
    for s1 in range(2):
        for s2 in range(2):
            result += sigmas_func(SIGMAS[s1], SIGMAS[s2], *sigmas_args)*terms[s1, s2]
    return result


@nb.njit
def cF2Function_norm(h, l):
    _terms = _Z2_E_terms_norm(h, l)*_Z2_F_terms(h, l)
    def _sigmas_func_12(s1, s2, *sigmas_args): return s1*s2
    return _Z2_sum_terms(_terms, _sigmas_func_12, 1, 1)


@nb.njit
def EFunction_norm(h, l, arg_s1, arg_s2):
    def _sigmas_func(s1, s2, arg_s1, arg_s2): return (s1 if arg_s1 else 1)*(s2 if arg_s2 else 1)
    return _Z2_sum_terms(_Z2_E_terms_norm(h, l), _sigmas_func, arg_s1, arg_s2)


# ======================================== MOMENTS ========================================
M_SQRT2 = np.sqrt(2)


def calc_Z2_norm(h, l):
    return M_SQRT2*(2*np.pi)**2/(h*np.sqrt(l))*cF2Function_norm(h, l)


def calc_m(h, l):
    if h == 0:
        return 0
    if l == 0:
        return Z1.calc_m_single(h)
    return -1/(2*h) - h/(2*l) + M_SQRT2*EFunction_norm(h, l, 1, 1)/(4*np.sqrt(l)*cF2Function_norm(h, l))


def calc_eta(h, l):
    if h == 0:
        return Z1.calc_m_single(l)
    if l == 0:
        return Z1.calc_m_single(h)**2
    return -1 - 1/(2*l) + h**2/(2*l**2) + M_SQRT2*EFunction_norm(h, l, 1, 0)/(2*np.sqrt(l)*cF2Function_norm(h, l)) - M_SQRT2*h*EFunction_norm(h, l, 1, 1)/(4*l**(3/2)*cF2Function_norm(h, l))


def calc_upsilon(m, eta, h, l):
    if h == 0 and l == 0:
        return 1/3
    if h == 0:
        return eta/(l*(eta+1))
    if l == 0:
        return m/h
    m_eta = -eta*(l + h**2)/(2*l*h) - m/l + (2*l*m*h + l + h**2)*EFunction_norm(h, l, 0, 0)/(2*l*h*EFunction_norm(h, l, 1, 1)) - 1/(2*l*h)
    return .5*(1-m_eta/m)


def calc_eta2(m, eta, h, l):
    if h == 0 and l == 0:
        return 1/3
    if h == 0:
        return 1-2*eta/l
    if l == 0:
        m_h = m/h
        return 6*m_h*m_h - 4*m_h + 1
    return eta*(-3*l + h**2)/(2*l**2) + m*h*(l + 1)/l**2 + (2*l**2 + 1)/(2*l**2) - (2*l*m*h + l + h**2)*EFunction_norm(h, l, 0, 0)/(2*l**2*EFunction_norm(h, l, 1, 1))


_, calc_Z2_norm = Z1.single_and_vectorize(calc_Z2_norm)
_, calc_m = Z1.single_and_vectorize(calc_m)
_, calc_eta = Z1.single_and_vectorize(calc_eta)
calc_upsilon_single, calc_upsilon = Z1.single_and_vectorize(calc_upsilon)
calc_eta2_single, calc_eta2 = Z1.single_and_vectorize(calc_eta2)


def calc_m_par_2(m, eta, h, l):
    if h == 0:
        m_h = (1+eta)/3
        return 1 - 2*m_h + 2*m_h*eta/(eta+1)
    if l == 0:
        return 1 - 2*m/h
    upsilon = calc_upsilon_single(m, eta, h, l)
    return 1 - 2*m/h + 2*l*m/h*upsilon


def calc_psi0(m, eta, h, l):
    upsilon = calc_upsilon_single(m, eta, h, l)
    eta2 = calc_eta2_single(m, eta, h, l)
    if h == 0:
        m_h = (1+eta)/3
    else:
        m_h = m/h
    return 1.5*m_h*upsilon + .25*(eta2-1)


_, calc_m_par_2 = Z1.single_and_vectorize(calc_m_par_2)
_, calc_psi0 = Z1.single_and_vectorize(calc_psi0)


# ======================================== FINDING ========================================
@nb.njit
def m_upsilon_approx(m, eta):
    m2 = m * m
    zeta = (eta - m2) / (1 - m2)
    zeta2 = zeta * zeta

    mu_p = Z1.mu_p_approx(m)
    m_upsilon = (1 - eta) / (1 - m2) * m * mu_p + 0.0467246 * m * (1 - m) * zeta * (1 - zeta) * (1 - 7.96819 * m - 0.939775 * zeta + 14.0242 * m * zeta + 20.9323 * m2 + 7.65054 * zeta2 - 14.5937 * m2 * zeta - 5.35767 * m * zeta2 - 13.542 * m * m2 - 3.72091 * zeta * zeta2) + (0.00349416 + (0.0191611 + (-0.0367299 + (-0.0684992 + (0.149972 - 0.067404 * m2) * m2) * m2) * m2) * m2) * m
    return m_upsilon


@nb.njit
def rho_approx(m, eta):
    m2 = m * m
    zeta = (eta - m2) / (1 - m2)
    rho = 1 - .5 * zeta - .514 * (1 - zeta) * zeta * (1 - (.232 + (1.681 - 1.466 * m) * m) * m + (-.123 * m - 1.041 + (0.858 - 0.356 * zeta) * zeta) * zeta)
    return rho


@nb.njit
def coeffs_approx(m, eta):
    mu_p = Z1.mu_p_approx(m)
    upsilon = m_upsilon_approx(m, eta) / m
    rho = rho_approx(m, eta)

    return rho * m / mu_p, (1 - rho) / upsilon


@nb.njit
def find_coeffs_single(m, eta, eps=1e-8):
    if eta < m * m:
        return np.nan, np.nan

    if m == 1 or eta == 1:
        return np.nan, np.nan
    
    if eta == m * m:
        return Z1.find_coeffs_single(m), 0

    if m == 0 and eta == 0:
        return 0, 0
    
    if m == 0:
        return 0, Z1.find_coeffs_single(eta)


    h, l = coeffs_approx(m, eta)

    delta_h, delta_l = 1, 1
    steps = 0
    while delta_h * delta_h + delta_l * delta_l > eps * eps and steps < 64:
        if h == 0 or l == 0:
            break
        if l < 0:
            break

        _l = 1 / l
        _h = 1 / h
        _sqrt_l = np.sqrt(_l)
        _exp_2l = np.exp(-2 * l)
        _exp_2h = np.exp(-2 * h)
        cF2     = _exp_2h * _exp_2h * F2((.5) * M_SQRT2 * _sqrt_l * (h - 2 * l)) - 2 * _exp_2h * _exp_2l * F2((.5) * M_SQRT2 * _sqrt_l * h) + F2((.5) * M_SQRT2 * _sqrt_l * (h + 2 * l))
        _cF2 = 1 / cF2
        E12_cF2 = _cF2 * (_exp_2h * _exp_2h - 2 * _exp_2h * _exp_2l + 1)
        E1_cF2 = _cF2 * (1 - _exp_2h * _exp_2h)
        E_cF2 = _cF2 * (_exp_2h * _exp_2h + 2 * _exp_2h * _exp_2l + 1)

        expr_m = (.25) * M_SQRT2 * E12_cF2 * _sqrt_l - .5 * _h - .5 * _l * h
        expr_m_diff_h = (.25) * M_SQRT2 * E12_cF2 * _l * _sqrt_l * h + (.5) * M_SQRT2 * E1_cF2 * _sqrt_l + (.5) * _h * _h + _l * (-.25 * E12_cF2 * E12_cF2 - 1.0 / 2.0)
        expr_m_diff_l = -.25 * E12_cF2 * E1_cF2 * _l - 0.125 * M_SQRT2 * E12_cF2 * _l * _l * _sqrt_l * h * h - 0.125 * M_SQRT2 * E12_cF2 * _l * _sqrt_l + _l * _l * ((0.125) * E12_cF2 * E12_cF2 * h + (.5) * h) + _sqrt_l * ((.25) * M_SQRT2 * E12_cF2 + (.25) * M_SQRT2 * E_cF2)

        expr_eta = -.25 * M_SQRT2 * E12_cF2 * _l * _sqrt_l * h + (.5) * M_SQRT2 * E1_cF2 * _sqrt_l + (.5) * _l * _l * h * h - .5 * _l - 1
        expr_eta_diff_h = -.5 * E12_cF2 * E1_cF2 * _l - .25 * M_SQRT2 * E12_cF2 * _l * _l * _sqrt_l * h * h - .25 * M_SQRT2 * E12_cF2 * _l * _sqrt_l + _l * _l * ((.25) * E12_cF2 * E12_cF2 * h + h) + _sqrt_l * ((.5) * M_SQRT2 * E12_cF2 + (.5) * M_SQRT2 * E_cF2)
        expr_eta_diff_l = (0.125) * M_SQRT2 * E12_cF2 * _l * _l * _l * _sqrt_l * h * h * h - .5 * E1_cF2 * E1_cF2 * _l + M_SQRT2 * E1_cF2 * _sqrt_l + _l * _l * _l * (-0.125 * E12_cF2 * E12_cF2 * h * h - h * h) + _l * _l * _sqrt_l * ((0.375) * M_SQRT2 * E12_cF2 * h - .25 * M_SQRT2 * E1_cF2 * h * h) + _l * _l * ((.5) * E12_cF2 * E1_cF2 * h + .5) + _l * _sqrt_l * (-.25 * M_SQRT2 * E12_cF2 * h - .25 * M_SQRT2 * E1_cF2 - .25 * M_SQRT2 * E_cF2 * h)

        D = expr_m_diff_h * expr_eta_diff_l - expr_m_diff_l * expr_eta_diff_h
        if D == 0:
            break
        _D = 1 / D
        delta_h = ((expr_m - m) * expr_eta_diff_l - (expr_eta - eta) * expr_m_diff_l) * _D
        delta_l = ((expr_eta - eta) * expr_m_diff_h - (expr_m - m) * expr_eta_diff_h) * _D

        h -= delta_h
        l -= delta_l
        steps += 1

    if delta_h * delta_h + delta_l * delta_l > eps * eps:
        return np.nan, np.nan

    return h, l


@nb.njit(parallel=True)
def find_coeffs_vectorize(m, eta, eps=1e-8):
    h = np.empty_like(m)
    l = np.empty_like(m)

    for i in nb.prange(m.size):
        result = find_coeffs_single(m.flat[i], eta.flat[i], eps)
        h.flat[i] = result[0]
        l.flat[i] = result[1]

    return h, l


def find_coeffs(m, eta, eps=1e-8):
    if isinstance(m, np.ndarray):
        return find_coeffs_vectorize(m, eta, eps)
    else:
        return find_coeffs_single(m, eta, eps)
