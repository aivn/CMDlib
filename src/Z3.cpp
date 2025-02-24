#include <math.h>

#include "../include/Z3.hpp"
#include "../include/special.hpp"

const int SIGMAS[2] = { -1, 1 };

// =======================================================================================

double cF3_norm_symmetrical(double h_i, double h_j, double l_ij, double l_ik) {
    double result = 0;
    for (int is1 = 0; is1 < 2; is1++)
        for (int is2 = 0; is2 < 2; is2++)
            for (int is3 = 0; is3 < 2; is3++) {
                int s1 = SIGMAS[is1];
                int s2 = SIGMAS[is2];
                int s3 = SIGMAS[is3];

                double arg_exp = (s3 - 1) * h_j + (s1 * s2 - 1) * l_ik + (s1 + s2) * (h_i + s3 * l_ij) - 2 * (h_i + l_ij);
                double arg1    = (h_i + s3 * l_ij + (s1 + s2) * h_i / h_j * l_ij) / sqrt(2 * h_i / h_j * l_ij);
                double arg2    = sqrt(h_i / h_j * l_ij / l_ik);
                double arg3    = (s1 + s2) / sqrt(2 * l_ik) * (l_ik - h_i / h_j * l_ij);
                result += s1 * s2 * s3 * exp(arg_exp) * F3(arg1, arg2, arg3);
            }

    return result;
}

double cF3_symmetrical(double h_i, double h_j, double l_ij, double l_ik) {
    return cF3_norm_symmetrical(h_i, h_j, l_ij, l_ik) * exp(h_j + l_ik + 2 * (h_i + l_ij));
}

double E_norm_symmetrical(double h_i, double h_j, double l_ij, double l_ik, bool arg_s1, bool arg_s2, bool arg_s3) {
    double result = 0;
    for (int is1 = 0; is1 < 2; is1++)
        for (int is2 = 0; is2 < 2; is2++)
            for (int is3 = 0; is3 < 2; is3++) {
                int s1 = SIGMAS[is1];
                int s2 = SIGMAS[is2];
                int s3 = SIGMAS[is3];

                int factor = 1;
                factor *= arg_s1 ? s1 : 1;
                factor *= arg_s2 ? s2 : 1;
                factor *= arg_s3 ? s3 : 1;

                double arg_exp = (s3 - 1) * h_j + (s1 * s2 - 1) * l_ik + (s1 + s2) * (h_i + s3 * l_ij) - 2 * (h_i + l_ij);
                result += factor * exp(arg_exp);
            }

    return result;
}

double E_symmetrical(double h_i, double h_j, double l_ij, double l_ik, bool arg_s1, bool arg_s2, bool arg_s3) {
    return E_norm_symmetrical(h_i, h_j, l_ij, l_ik, arg_s1, arg_s2, arg_s3) * exp(h_j + l_ik + 2 * (h_i + l_ij));
}

double FTilde_norm_symmetrical(double h_i, double h_j, double l_ij, double l_ik, bool arg_s1, bool arg_s2, bool arg_s3) {
    double result = 0;
    for (int is1 = 0; is1 < 2; is1++)
        for (int is2 = 0; is2 < 2; is2++)
            for (int is3 = 0; is3 < 2; is3++) {
                int s1 = SIGMAS[is1];
                int s2 = SIGMAS[is2];
                int s3 = SIGMAS[is3];

                int factor = 1;
                factor *= arg_s1 ? s1 : 1;
                factor *= arg_s2 ? s2 : 1;
                factor *= arg_s3 ? s3 : 1;

                double arg_exp    = (s3 - 1) * h_j + (s1 * s2 - 1) * l_ik + (s1 + s2) * (h_i + s3 * l_ij) - 2 * (h_i + l_ij);
                double arg_dawson = (h_i + s3 * l_ij + (s1 + s2) * l_ik) / sqrt(2 * l_ik);
                result += factor * exp(arg_exp) * dawson(arg_dawson);
            }

    return result;
}

double FTilde_symmetrical(double h_i, double h_j, double l_ij, double l_ik, bool arg_s1, bool arg_s2, bool arg_s3) {
    return FTilde_norm_symmetrical(h_i, h_j, l_ij, l_ik, arg_s1, arg_s2, arg_s3) * exp(h_j + l_ik + 2 * (h_i + l_ij));
}

double F_norm_symmetrical(double h_i, double h_j, double l_ij, double l_ik, bool arg_s1, bool arg_s2, bool arg_s3) {
    double result = 0;
    for (int is1 = 0; is1 < 2; is1++)
        for (int is2 = 0; is2 < 2; is2++)
            for (int is3 = 0; is3 < 2; is3++) {
                int s1 = SIGMAS[is1];
                int s2 = SIGMAS[is2];
                int s3 = SIGMAS[is3];

                int factor = 1;
                factor *= arg_s1 ? s1 : 1;
                factor *= arg_s2 ? s2 : 1;
                factor *= arg_s3 ? s3 : 1;

                double arg_exp    = (s3 - 1) * h_j + (s1 * s2 - 1) * l_ik + (s1 + s2) * (h_i + s3 * l_ij) - 2 * (h_i + l_ij);
                double arg_dawson = (h_i + s3 * l_ij + (s1 + s2) * h_i / h_j * l_ij) / sqrt(2 * h_i / h_j * l_ij);
                result += factor * exp(arg_exp) * dawson(arg_dawson);
            }

    return result;
}

double F_symmetrical(double h_i, double h_j, double l_ij, double l_ik, bool arg_s1, bool arg_s2, bool arg_s3) {
    return F_norm_symmetrical(h_i, h_j, l_ij, l_ik, arg_s1, arg_s2, arg_s3) * exp(h_j + l_ik + 2 * (h_i + l_ij));
}

// =======================================================================================
