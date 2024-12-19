#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Z2.hpp"

inline double Z(double p){ return 4*M_PI*(fabs(p)>1e-4? sinh(p)/p: 1+p*p/6); }
inline double L(double p){ return fabs(p)>1e-4? 1/tanh(p)-1/p: p/3; }
inline double dLdp(double p){ double shp = sinh(p); return fabs(p)>1e-6? 1/(p*p)-1/(shp*shp) : 1/3.-p*p/15; }
// const float MAX_P = logf(HUGE_VALF)-1;
inline double invL(double M){
	if(fabsf(M)<1e-6) return M*3;
	double p = 1.f;
	for(int i=0; i<10; i++){
		double Lp = L(p);
		if(p>400 || fabs(Lp-M)<1e-6) break;
		p = p-(Lp-M)/dLdp(p);
		if(fabs(p)<1e-6) return p;
	}
	return p;
}

int main(int argc, const char **argv){
	if(argc!=3){ printf("usage: ./calc-M-eta-to-rho-p-lambda M eta ==> rho p lambda Z2 S1 S2\n"); exit(1); }
	double M = atof(argv[1]), eta = atof(argv[2]), p = invL(M);
	Z2 z2; z2.calc(M, eta); double z = z2.calc_Z2(), S2 = log(z)/2 - z2.h*M -z2.lambda*eta/2;
	printf("%f %f %f %f %f %f\n", (p>1e-3? z2.h/p: 1), p, (eta>M*M?z2.lambda:0), z, log(Z(p))-p*M, S2);
	
	return 0;
}
