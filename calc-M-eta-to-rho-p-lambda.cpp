#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "include/Z2.hpp"

inline double Z(double p){ return 4*M_PI*(fabs(p)>1e-4? sinh(p)/p: 1+p*p/6); }

double calc_mu_p(double M){
	double M2 = M*M, tmp = 0.775383f + (0.109185f + (0.289114f + (-0.871214f + (1.85968f + (-1.71306f + 0.550912f*M2 )*M2 )*M2 )*M2 )*M2 )*M2;
	return (1-M2*tmp*tmp)/3.;
}

int main(int argc, const char **argv){
	if(argc!=3){ printf("usage: ./calc-M-eta-to-rho-p-lambda M eta ==> rho p lambda Z2 S1 S2 eta2\n"); exit(1); }
	double M = atof(argv[1]), eta = atof(argv[2]), p = invL(M), zeta = (eta-M*M)/(1-M*M);
	Z2_symmetrical z2;  z2.m = M; z2.eta = eta;
	double rho = 1-zeta/2-0.514*(1-zeta)*zeta*( 1 - 0.232*M + 1.681*M*M - 1.466*M*M*M - 0.123*M*zeta
												- 1.041*zeta + 0.858*zeta*zeta - 0.356*zeta*zeta*zeta);
	
	z2.calc_from_moments(p*rho, (1-rho)/(3*calc_mu_p(M)*calc_mu_p(zeta)/(1+zeta))); double z = z2.Z2_norm, S2 = log(z)/2 - z2.h*M -z2.l*eta/2;
	printf("%f %f %f %f %f %f %f\n", (p>1e-3? z2.h/p: 1), p, (eta>M*M?z2.l:0), z, log(Z(p))-p*M, S2, z2.eta2);
	
	return 0;
}
