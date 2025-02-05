#include <cassert>
#include <iostream>

#include <math.h>

#include "../include/Z2.hpp"

using namespace Z2_symmetrical;

int main() {
    Z2 z2 = Z2();

    // z2.h = 1;
    // z2.l = 2;
    // z2.calc_from_coeffs();
    // assert(fabs(z2.m - 0.4386289825963412) < 1e-6);
    // assert(fabs(z2.eta - 0.5836393087685698) < 1e-6);
    // assert(fabs3(z2.upsilon - 0.16575974898177726) < 1e-6);
    // assert(fabs(z2.eta2 - 0.489067721282732) < 1e-6);
    // assert(fabs(z2.m_par_2 - 0.41357015501252503) < 1e-6);
    // assert(fabs(z2.psi0 + 0.01867252460236421) < 1e-6);

    // z2.calc_from_moments(.5, .5);
    // assert(fabs(z2.h - 1) < 1e-7);
    // assert(fabs(z2.l - 2) < 1e-7);

    // z2.m   = 0.5;
    // z2.eta = 0.1;
    // z2.calc_from_moments(.5, .5);
    // assert(std::isnan(z2.h));
    // assert(std::isnan(z2.l));
    // assert(std::isnan(z2.upsilon));

    double min  = 0;
    double max  = 1;
    double N    = 100;
    double step = (max - min) / N;

    // 0.0107152 0.0912011 0.00368028 0.0304131 0.430377
    std::cout << "#:h l m eta zeta steps upsilon eta2 m_par_2 psi0 Z2_norm" << '\n';
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            double m = min + j * step, zeta = min + i * step;
            double eta = zeta + (1 - zeta) * m * m;

            double m2    = m * m;
            double zeta2 = zeta * zeta;

            double mu_p_arg = m * ((0.775383) + ((0.109185) + ((0.289114) + ((-0.871214) + ((1.85968) + ((-1.71306) + (0.550912) * m2) * m2) * m2) * m2) * m2) * m2);
            double mu_p     = (1 - mu_p_arg * mu_p_arg) / 3;

            double upsilonM = (1 - eta) / (1 - m2) * m * mu_p + 0.0467246 * m * (1 - m) * zeta * (1 - zeta) * (1 - 7.96819 * m - 0.939775 * zeta + 14.0242 * m * zeta + 20.9323 * m2 + 7.65054 * zeta2 - 14.5937 * m2 * zeta - 5.35767 * m * zeta2 - 13.542 * m * m2 - 3.72091 * zeta * zeta2) + (0.00349416 + (0.0191611 + (-0.0367299 + (-0.0684992 + (0.149972 - 0.067404 * m2) * m2) * m2) * m2) * m2) * m;
            double upsilon  = upsilonM / m;
            double rho      = 1 - .5 * zeta - .514 * (1 - zeta) * zeta * (1 - (.232 + (1.681 - 1.466 * m) * m) * m + (-.123 * m - 1.041 + (0.858 - 0.356 * zeta) * zeta) * zeta);

            z2.m   = m;
            z2.eta = eta;
            int steps = z2.calc_from_moments(rho * m / mu_p, (1 - rho) / upsilon);
            std::cout << z2.h << ' ' << z2.l << ' ';
            std::cout << z2.m << ' ' << z2.eta << ' ' << zeta << ' ';
            std::cout << steps << ' ';
            std::cout << z2.upsilon << ' ' << z2.eta2 << ' ' << z2.m_par_2 << ' ' << z2.psi0 << ' ' << z2.Z2_norm;
            std::cout << '\n';
        }
        std::cout << '\n';
    }
    return 0;
}
