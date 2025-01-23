#ifndef SPEC_FUNC_HPP
#define SPEC_FUNC_HPP

#include "aitken.hpp"
#include <complex>

class SpecFunc{
public:
    template <typename T1, typename T2> T1 Infeld_func(T1 x, T2 nu, T1 eps){
        T1 res = aitken([nu, x](T1 xi){return exp(x*cos(xi))*pow(sin(xi), 2*nu);}, 0., M_PI, eps);
        return res*pow(2, -nu)*pow(x, nu)/sqrt(M_PI)/tgamma(nu+0.5);
    }
};

#endif //SPEC_FUNC_HPP