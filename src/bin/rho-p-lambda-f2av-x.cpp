#include <stdlib.h>
#include <math.h>
#include <iostream>

inline double lang(double p){ return fabs(p>1e-6)? 1/tanh(p)-1/p: p/3; }
inline double lang2(double p){ return fabs(p>1e-6)? 2/(p*p) - 2/(p*tanh(p)) +1: 1./3; }
inline double pow2(double x){ return x*x; }

double calc_p(double M){
	if(fabs(M<1e-6)) return M*3;
	double p = 1.;
	for(int i=0; i<10; i++){
		double sp = sinh(p);
		p = p-(lang(p)-M)/(1/(p*p) -1/(sp*sp));
		if(fabs(p)<1e-6) return p;
		if(p>400 || fabs(lang(p)-M)<1e-8) break;
	}
	return p;
}

inline double Zfunc(double p){ return 4*M_PI*(p>1e-6? sinh(p)/p: 1.+p*p/6); }
inline double dZdp(double p){ return 4*M_PI*(fabs(p)>1e-4? (p*cosh(p)-sinh(p))/(p*p): p/3); }


int main(int argc, const char **argv){
	if(argc<2 || argc>4){ printf("usage: ./rho-p-lambda-f2av-x x_sz [M_sz [eta_sz]] ==> #:M eta M2 Meta  rho p lambda_ Upsilon zeta  Z\n"); return 1; }
	int x_sz = atof(argv[1]), M_sz = 50, eta_sz = 50;
	if(argc>2) M_sz = atoi(argv[2]);
	if(argc>3) eta_sz = atoi(argv[3]);	
	double dx = 2./x_sz, dM = 1./M_sz, deta = 1./eta_sz;
	std::cout<<"#:M eta M2 Meta  rho p lmbda Upsilon zeta  Z2\n";  std::cout.precision(17);

	for(int iM=0; iM<M_sz; iM++){
		double M = iM*dM,  rhop = calc_p(M);
		for(int ieta=0; ieta<eta_sz; ieta++){
			double eta = ieta*deta, v = calc_p(eta);
			if(iM==0){		
				//         M       eta       M2         Meta      rho             p      lmbda     Upsilon                          zeta      Z
				std::cout<<0<<' '<<eta<<' '<<1./3<<' '<<0<<"   "<<1/(1+eta)<<' '<<0<<' '<<v<<' '<<(v==0?1./3:eta/(v*(1+eta)))<<' '<<eta<<' '<<4*M_PI*Zfunc(v)<<std::endl;
				continue;
			}
			if(ieta==0){  // приближение среднего поля
				double  M2 = 1-2*M/rhop;  
				//         M       eta       M2       Meta         rho     p          lmbda   Upsilon    zeta
				std::cout<<M<<' '<<M*M<<' '<<M2<<' '<<M2*M<<"   "<<1<<' '<<rhop<<' '<<0<<' '<<M/rhop<<' '<<0<<' '<<pow2(Zfunc(rhop))<<std::endl;
				continue;
			}
			M = eta = 0; double Z = 0, Meta = 0, M2 = 0, A = rhop*rhop+v*v, B = 2*rhop*v;
			for(int ix=0; ix<x_sz; ix++){
				double x = (ix+.5)*dx - 1;
				double e = exp(rhop*x);
				double arg = sqrt(A+B*x);
				double z = Zfunc(arg);
				Z += e*z;
				M += x*e*z;
				M2 += x*x*e*z;
				double de = (v+rhop*x)*e*dZdp(arg)/arg;
				eta += de;
				Meta += x*de;
			}

			Z *= 2*M_PI*dx;
			M *= 2*M_PI*dx/Z;
			M2 *= 2*M_PI*dx/Z;
			eta *= 2*M_PI*dx/Z;
			Meta *= 2*M_PI*dx/Z;
			double p = calc_p(M);

			//         M       eta       M2       Meta         rho          p       lmbda   Upsilon            zeta
			std::cout<<M<<' '<<eta<<' '<<M2<<' '<<Meta<<"   "<<rhop/p<<' '<<p<<' '<<v<<' '<<(1-rhop/p)/v<<' '<<(eta-M*M)/(1-M*M)<<' '<<Z<<std::endl;
		}
		M = lang(2*rhop);  // eta = 1
		//         M       eta     M2                   Meta      rho          p     lmbda   Upsilon  zeta
		std::cout<<M<<' '<<1<<' '<<1-2*M/(2*rhop)<<' '<<M<<"   "<<.5<<' '<<(2*rhop)<<" inf "<<0<<' '<<1<<' '<<"inf"<<'\n'<<std::endl;
	}
	
	for(int ieta=0; ieta<eta_sz; ieta++){ // M=1
		double eta = ieta*deta, v = calc_p(eta);
		//         M       eta     M2      Meta       rho  p    lmbda     Upsilon zeta
		std::cout<<1<<' '<<1<<' '<<1<<' '<<1<<"   "<<.5<<" inf "<<v<<' '<<0<<' '<<eta<<' '<<"inf"<<std::endl;
	}
	//         M       eta     M2      Meta       rho  p   lmbda  Upsilon zeta
	std::cout<<1<<' '<<1<<' '<<1<<' '<<1<<"   "<<.5<<" inf inf "<<0<<' '<<1<<' '<<"inf"<<std::endl;
	
	return 0;
}
