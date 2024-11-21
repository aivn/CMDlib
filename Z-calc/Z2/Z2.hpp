#include <math.h>

#include "dawson/gsl_sf_dawson.h"

const double ACCURACY = 1e-8;

struct Z2 {
private:
    double E12_cF2, E1_cF2, E_cF2;

public:
    double h, lambda;
    int steps;

    void calc(double m, double eta) {
        double m2    = m * m;
        double zeta  = (eta - m2) / (1 - m2);
        double zeta2 = zeta * zeta;

        double mu_p_arg = m * ((0.775383) + ((0.109185) + ((0.289114) + ((-0.871214) + ((1.85968) + ((-1.71306) + (0.550912) * m2) * m2) * m2) * m2) * m2) * m2);
        double mu_p     = (1 - mu_p_arg * mu_p_arg) / 3;

        double upsilonM = (1 - eta) / (1 - m2) * m * mu_p + 0.0467246 * m * (1 - m) * zeta * (1 - zeta) * (1 - 7.96819 * m - 0.939775 * zeta + 14.0242 * m * zeta + 20.9323 * m2 + 7.65054 * zeta2 - 14.5937 * m2 * zeta - 5.35767 * m * zeta2 - 13.542 * m * m2 - 3.72091 * zeta * zeta2) + (0.00349416 + (0.0191611 + (-0.0367299 + (-0.0684992 + (0.149972 - 0.067404 * m2) * m2) * m2) * m2) * m2) * m;
        double upsilon  = upsilonM / m;
        double rho      = 1 - .5 * zeta - .514 * (1 - zeta) * zeta * (1 - (.232 + (1.681 - 1.466 * m) * m) * m + (-.123 * m - 1.041 + (0.858 - 0.356 * zeta) * zeta) * zeta);

        h      = rho * m / mu_p;
        lambda = (1 - rho) / upsilon;

        double delta_h = 1, delta_l = 1;
        steps = 0.;
        while (delta_h * delta_h + delta_l * delta_l > ACCURACY * ACCURACY && steps < 64) {
            double _l      = 1 / lambda;
            double _h      = 1 / h;
            double _sqrt_l = sqrt(_l);
            double _exp_2l = exp(-2 * lambda);
            double _exp_2h = exp(-2 * h);
            double cF2     = _exp_2h * _exp_2h * gsl_sf_dawson((.5) * M_SQRT2 * _sqrt_l * (h - 2 * lambda)) - 2 * _exp_2h * _exp_2l * gsl_sf_dawson((.5) * M_SQRT2 * _sqrt_l * h) + gsl_sf_dawson((.5) * M_SQRT2 * _sqrt_l * (h + 2 * lambda));
            double _cF2    = 1 / cF2;
            E12_cF2        = _cF2 * (_exp_2h * _exp_2h - 2 * _exp_2h * _exp_2l + 1);
            E1_cF2         = _cF2 * (1 - _exp_2h * _exp_2h);
            E_cF2          = _cF2 * (_exp_2h * _exp_2h + 2 * _exp_2h * _exp_2l + 1);

            double G1        = (.25) * M_SQRT2 * E12_cF2 * _sqrt_l - .5 * _h - .5 * _l * h - m;
            double G1_diff_h = (.25) * M_SQRT2 * E12_cF2 * _l * _sqrt_l * h + (.5) * M_SQRT2 * E1_cF2 * _sqrt_l + (.5) * _h * _h + _l * (-.25 * E12_cF2 * E12_cF2 - 1.0 / 2.0);
            double G1_diff_l = -.25 * E12_cF2 * E1_cF2 * _l - 0.125 * M_SQRT2 * E12_cF2 * _l * _l * _sqrt_l * h * h - 0.125 * M_SQRT2 * E12_cF2 * _l * _sqrt_l + _l * _l * ((0.125) * E12_cF2 * E12_cF2 * h + (.5) * h) + _sqrt_l * ((.25) * M_SQRT2 * E12_cF2 + (.25) * M_SQRT2 * E_cF2);

            double G2        = -.25 * M_SQRT2 * E12_cF2 * _l * _sqrt_l * h + (.5) * M_SQRT2 * E1_cF2 * _sqrt_l + (.5) * _l * _l * h * h - .5 * _l - eta - 1;
            double G2_diff_h = -.5 * E12_cF2 * E1_cF2 * _l - .25 * M_SQRT2 * E12_cF2 * _l * _l * _sqrt_l * h * h - .25 * M_SQRT2 * E12_cF2 * _l * _sqrt_l + _l * _l * ((.25) * E12_cF2 * E12_cF2 * h + h) + _sqrt_l * ((.5) * M_SQRT2 * E12_cF2 + (.5) * M_SQRT2 * E_cF2);
            double G2_diff_l = (0.125) * M_SQRT2 * E12_cF2 * _l * _l * _l * _sqrt_l * h * h * h - .5 * E1_cF2 * E1_cF2 * _l + M_SQRT2 * E1_cF2 * _sqrt_l + _l * _l * _l * (-0.125 * E12_cF2 * E12_cF2 * h * h - h * h) + _l * _l * _sqrt_l * ((0.375) * M_SQRT2 * E12_cF2 * h - .25 * M_SQRT2 * E1_cF2 * h * h) + _l * _l * ((.5) * E12_cF2 * E1_cF2 * h + .5) + _l * _sqrt_l * (-.25 * M_SQRT2 * E12_cF2 * h - .25 * M_SQRT2 * E1_cF2 - .25 * M_SQRT2 * E_cF2 * h);

            double _D = 1 / (G1_diff_h * G2_diff_l - G1_diff_l * G2_diff_h);
            delta_h   = (G1 * G2_diff_l - G2 * G1_diff_l) * _D;
            delta_l   = (G2 * G1_diff_h - G1 * G2_diff_h) * _D;

            h -= delta_h;
            lambda -= delta_l;
            steps += 1;
        }

        if (h != h || lambda != lambda) {
            h      = std::nan("1");
            lambda = std::nan("1");
            steps  = std::nan("1");
        }
    }
};
