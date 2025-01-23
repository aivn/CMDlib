#ifndef COP_HPP
#define COP_HPP

#include <vector>
#include "polynomials.hpp"
#include <iostream>
#include <complex>

std::complex<double> iu(0, 1);

int factorial(int n){
    if(n==0){
        return 1;
    }else{
        return n*factorial(n-1);
    }
}

class COP{
    std::vector<Polynomial<double>> Legendre_Pol;
    std::vector<std::vector<Polynomial<double>>> Legendre_Pol_derivatives;
    int n;
public:

    int get_n() const { return n; }

    void init_Legendre_Pol(int n_){
        n = n_;
        Polynomial<double> p1(0); p1[0] = 1; Polynomial<double> p2(1); p2[1] = 1;
        Legendre_Pol.push_back(p1); Legendre_Pol.push_back(p2);

        Polynomial<double> x_pol(1); x_pol[1] = 1;
        for(int i = 2; i < n_+1; i++){
            Polynomial<double> p(i);
            p = double(2*i-1)/double(i)*x_pol*Legendre_Pol[i-1]-double(i-1)/double(i)*Legendre_Pol[i-2];
            Legendre_Pol.push_back(p); 
        }
    }

    Polynomial<double> get_Legendre_pol(int i) const{
        return Legendre_Pol[i];
    }

    void calc_Legendre_Pol_derivatives(){
        Legendre_Pol_derivatives.resize(n+1);
        for(int l = 0; l < n+1; l++){
            Polynomial<double> p = Legendre_Pol[l];
            for(int m = 0; m <= l; m++){
                Legendre_Pol_derivatives[l].push_back(p.diff(m));
            }
        }
    }

    std::complex<double> Y(int l, int m, double teta, double fi) const{
        int m_abs = std::abs(m);
        if(m < 0){
            return pow(-1, m_abs)*std::conj(Y(l, m_abs, teta, fi));
        }else{
            return pow(-1, m)*sqrt(double(2*l+1)*double(factorial(l-m))/(4*M_PI*double(factorial(l+m))))*exp(fi*iu*double(m))*Legendre_Pol_derivatives[l][m].calc(cos(teta))*pow(1-cos(teta)*cos(teta), double(m)/2);
        }
    }
};

#endif //COP_HPP