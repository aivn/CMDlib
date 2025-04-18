#include <math.h>

#include "../include/Z2.hpp"
#include "../include/special.hpp"

// #include <iostream>

const int SIGMAS[2] = { -1, 1 };

// =======================================================================================

double cF2_norm_symmetrical(double h, double l) {
    double result = 0;
    for (int is1 = 0; is1 < 2; is1++)
        for (int is2 = 0; is2 < 2; is2++) {
            int s1 = SIGMAS[is1];
            int s2 = SIGMAS[is2];
            result += s1 * s2 * exp((s1 + s2 - 2) * h + (s1 * s2 - 1) * l) * dawson((h + (s1 + s2) * l) / sqrt(2 * l));
        }

    return result;
}

double cF2_symmetrical(double h, double l) {
    return cF2_norm_symmetrical(h, l) * exp(2 * h + l);
}

double E_norm_symmetrical(double h, double l, bool arg_s1, bool arg_s2) {
    double result = 0;
    for (int is1 = 0; is1 < 2; is1++)
        for (int is2 = 0; is2 < 2; is2++) {
            int s1 = SIGMAS[is1];
            int s2 = SIGMAS[is2];

            int factor = 1;
            factor *= arg_s1 ? s1 : 1;
            factor *= arg_s2 ? s2 : 1;

            result += factor * exp((s1 + s2 - 2) * h + (s1 * s2 - 1) * l);
        }

    return result;
}

double E_symmetrical(double h, double l, bool arg_s1, bool arg_s2) {
    return E_norm_symmetrical(h, l, arg_s1, arg_s2) * exp(2 * h + l);
}

// =======================================================================================

void Z2_symmetrical::calc_base() {
    // _h = fabs(h) > 1e-4 ? 1 / h : NAN;
    // _l = l > 1e-4 ? 1 / l : NAN;
    _h = 1 / h;
    _l = 1 / l;

    sqrt_l  = sqrt(l);
    _sqrt_l = 1 / sqrt_l;
    // _sqrt_l = sqrt_l > 1e-4 ? 1 / sqrt_l : NAN;

    _exp_2l = exp(-2 * l);
    _exp_2h = exp(-2 * h);

    m_h = NAN;
    // if (!std::isnan(_h) || !std::isnan(_l)) {
    cF2_norm = _exp_2h * _exp_2h * dawson((.5) * M_SQRT2 * _sqrt_l * (h - 2 * l)) - 2 * _exp_2h * _exp_2l * dawson((.5) * M_SQRT2 * _sqrt_l * h) + dawson((.5) * M_SQRT2 * _sqrt_l * (h + 2 * l));

    E12_norm = _exp_2h * _exp_2h - 2 * _exp_2h * _exp_2l + 1;
    E1_norm  = 1 - _exp_2h * _exp_2h;
    E_norm   = _exp_2h * _exp_2h + 2 * _exp_2h * _exp_2l + 1;
    // } else sqrt_l = _sqrt_l = cF2_norm = E12_norm = E1_norm = E_norm = NAN;

    // m_llb    = fabs(h) > 1e-4 ? (1 + _exp_2h) / (1 - _exp_2h) - _h : h / 3;
    // m_h_llb  = fabs(h) > 1e-4 ? m_llb * _h : 1. / 3.;
    // eta_prmg = l > 1e-4 ? (1 + _exp_2l) / (1 - _exp_2l) - _l : l / 3;
}

void Z2_symmetrical::calc_m() {
    // if (fabs(h) > 1e-2 && l > 1e-1) {
    m   = -.25 * (-M_SQRT2 * E12_norm * h * l + 2 * cF2_norm * h * h * sqrt_l + 2 * cF2_norm * l * sqrt_l) / (cF2_norm * h * l * sqrt_l);
    m_h = m * _h;
    // } else if (l <= 1e-1) {
    //     m   = m_llb + l * (m_llb - m_llb * m_llb * m_llb - 2 * m_llb * m_h_llb);
    //     m_h = m_h_llb + l * (m_h_llb - m_h_llb * m_llb * m_llb - 2 * m_h_llb * m_h_llb);
    // } else {  // fabs(h) <= 1e-2
    //     m_h = (1 + eta_prmg) / 3;
    //     m   = h * m_h;
    // }
}

void Z2_symmetrical::calc_eta() {
    // if (fabs(h) > 1e-2 && l > 1e-1)
    eta = .25 * (-M_SQRT2 * E12_norm * h * l * l + 2 * M_SQRT2 * E1_norm * l * l * l + 2 * cF2_norm * h * h * l * sqrt_l - 4 * cF2_norm * l * l * l * sqrt_l - 2 * cF2_norm * l * l * sqrt_l) / (cF2_norm * l * l * l * sqrt_l);
    // else if (l <= 1e-1) eta = m_llb * m_llb + l * (1 - m_llb * m_llb * m_llb * m_llb + 6 * m_h_llb * m_h_llb - 4 * m_h_llb);
    // else eta = eta_prmg;  // fabs(h) <= 1e-2
}

void Z2_symmetrical::calc_upsilon() {
    // if (l <= 1e-2 || (fabs(h) <= 1e-1 && l < .1)) {
    //     double tmp = fabs(h) > 1e-4 ? -9 * m_h_llb * _h * _h + 6 * _h * _h - _h / m_llb : 0;
    //     upsilon    = m_h_llb + l * (-m_llb * m_llb * m_h_llb + m_llb * m_llb / 2 - 4 * m_h_llb * m_h_llb + 3 * m_h_llb - .5 + tmp);
    //     m_eta      = m * (1 - 2 * upsilon);
    // } else if (fabs(h) <= 1e-2 || (fabs(h) <= 1e-1 && l < .5)) {
    //     upsilon = eta * _l / (eta + 1);
    //     m_eta   = m * (1 - 2 * upsilon);
    // } else {
    m_eta   = -.5 * (E12_norm * eta * (h * h) + E12_norm * eta * l + 2 * E12_norm * h * m + E12_norm - E_norm * (h * h) - 2 * E_norm * h * l * m - E_norm * l) / (E12_norm * h * l);
    upsilon = .5 * (1 - m_eta / m);
    // }
}

void Z2_symmetrical::calc_eta2() {
    // if (l > 1e-2)
    eta2 = 1 + 2 * _l * (h * m * upsilon - eta);
    // else if (l <= 1e-2) {
    //     double tmp = fabs(h) > 1e-1 ? 90 * m_h_llb * m_h_llb * _h * _h - 60 * m_h_llb * _h * _h + 10 * _h * _h : 0;
    //     eta2       = 6 * m_h_llb * m_h_llb - 4 * m_h_llb + 1 + l * (-6 * m_llb * m_llb * m_h_llb * m_h_llb + 4 * m_llb * m_llb * m_h_llb + 12 * m_h_llb * m_h_llb - 4 * m_h_llb + tmp);
    // }
}

void Z2_symmetrical::calc_mh2() {
    mh2 = 1 - 2 * m_h + 2 * l * m_h * upsilon;
}

void Z2_symmetrical::calc_psi0() {
    psi0 = 1.5 * m_h * upsilon + .25 * (eta2 - 1);
}

void Z2_symmetrical::calc_Z2_norm() {
    // if (fabs(h) > 1e-2 && l > 1e-2)
    Z2_norm = M_SQRT2 * (2 * M_PI) * (2 * M_PI) * cF2_norm / (h * sqrt_l);
    // else if (l <= 1e-2) {
    //     double Z = 4 * M_PI * (fabs(h) > 1e-4 ? .5 * (1 - _exp_2h) * _h : exp(-h) * (1 + h * h / 6));
    //     Z2_norm  = Z * Z * (1 + l * eta);
    // } else {  // fabs(h) <= 1e-2
    //     double Z = (4 * M_PI) * (4 * M_PI) * (l > 1e-4 ? .5 * (1 - _exp_2l) * _l : exp(-l) * (1 + l * l / 6));
    //     Z2_norm  = Z * (1 + h * m);
    // }
}

// =======================================================================================

void Z2_symmetrical::calc_from_coeffs() {
    calc_base();
    calc_m();
    calc_eta();
    calc_upsilon();
    calc_eta2();
    calc_mh2();
    calc_psi0();
    calc_Z2_norm();
    // if (eta < m * m || eta2 > 1) m = eta = upsilon = eta2 = mh2 = psi0 = NAN;
}

int Z2_symmetrical::calc_from_moments(double h0, double l0, double eps, int max_steps) {
    int steps = 0;

    if (eta < m * m) h = l = NAN;
    else if (m >= 1 || eta >= 1) h = l = NAN;
    else if (eta == m * m) {
        h = invL(m);
        l = 0;
    } else if (m == 0 && eta == 0) h = l = 0;
    else if (m == 0) {
        h = 0;
        l = invL(eta);
    } else {
        h = h0;
        l = l0;

        double m_ref = m, eta_ref = eta;

        double delta_h = 1, delta_l = 1;
        while (delta_h * delta_h + delta_l * delta_l > eps * eps && steps < max_steps) {
            if (h == 0 || l == 0) break;
            if (l < 0) break;

            // std::cout << h << ' ' << l << '\n';

            calc_base();

            calc_m();
            calc_eta();
            calc_upsilon();
            calc_eta2();

            double expr_m_diff_h = eta - 2 * m * m - 2 * m_h + 1;
            double expr_m_diff_l = m_eta - m * eta;

            double expr_eta_diff_h = 2 * expr_m_diff_l;
            double expr_eta_diff_l = eta2 - eta * eta;

            double D = expr_m_diff_h * expr_eta_diff_l - expr_m_diff_l * expr_eta_diff_h;
            if (D == 0) break;
            double _D = 1 / D;
            delta_h   = ((m - m_ref) * expr_eta_diff_l - (eta - eta_ref) * expr_m_diff_l) * _D;
            delta_l   = ((eta - eta_ref) * expr_m_diff_h - (m - m_ref) * expr_eta_diff_h) * _D;

            h -= delta_h;
            l -= delta_l;
            steps += 1;
        }

        if (h != h || l != l) h = l = NAN;
        if (delta_h * delta_h + delta_l * delta_l > eps * eps) h = l = NAN;
    }

    calc_from_coeffs();

    return steps;
}

// =======================================================================================

double cF2_norm_not_symmetrical(double hi, double hj, double l) {
    double result = 0;
    for (int is1 = 0; is1 < 2; is1++)
        for (int is2 = 0; is2 < 2; is2++) {
            int s1 = SIGMAS[is1];
            int s2 = SIGMAS[is2];
            result += s1 * s2 * exp((s1 - 1) * hi + (s2 - 1) * hj + (s1 * s2 - 1) * l) * dawson((hj + s1 * l + s2 * l * hj / hi) / sqrt(2 * l * hj / hi));
        }

    return result;
}

double cF2_not_symmetrical(double hi, double hj, double l) {
    return cF2_norm_not_symmetrical(hi, hj, l) * exp(hi + hj + l);
}

double E_norm_not_symmetrical(double hi, double hj, double l, bool arg_s1, bool arg_s2) {
    double result = 0;
    for (int is1 = 0; is1 < 2; is1++)
        for (int is2 = 0; is2 < 2; is2++) {
            int s1 = SIGMAS[is1];
            int s2 = SIGMAS[is2];

            int factor = 1;
            factor *= arg_s1 ? s1 : 1;
            factor *= arg_s2 ? s2 : 1;

            result += factor * exp((s1 - 1) * hi + (s2 - 1) * hj + (s1 * s2 - 1) * l);
        }

    return result;
}

double E_not_symmetrical(double hi, double hj, double l, bool arg_s1, bool arg_s2) {
    return E_norm_not_symmetrical(hi, hj, l, arg_s1, arg_s2) * exp(hi + hj + l);
}
