#include <omp.h>
#include <stdlib.h>
#include <aiwlib/sphere>
using namespace aiw;

inline double lang(double p){ return fabs(p>1e-6)? 1/tanh(p)-1/p: p/3; }
inline double lang2(double p){ return fabs(p>1e-6)? 2/(p*p) - 2/(p*tanh(p)) +1: 1./3; }
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

int main(int argc, const char **argv){
	double t0 = omp_get_wtime();
	if(argc<2 || argc>4){ printf("usage: ./rho-p-lambda-f2av rank  [M_sz [eta_sz]] ==> #:M eta eta2 M2 Meta  rho p lambda_ Upsilon zeta\n"); return 1; }
	int rank = atoi(argv[1]), M_sz = (argc>=3? atoi(argv[2]): 95), lambda_sz = (argc==4? atoi(argv[3]): 40); 
	sph_init_table(rank);	int sz = sph_vertex_num(rank);

	double dM = .95/M_sz, dv0 = 4./lambda_sz, dv1 = 1+dv0;
	
	std::cout<<"#:M eta eta2 M2 Meta  rho p lmbda Upsilon zeta\n";  std::cout.precision(17);

	for(double v=0; v<=20;){ // M=0
		double eta = lang(v);
		//         M       eta       eta2           M2         Meta      rho             p     lmbda     Upsilon
		std::cout<<0<<' '<<eta<<' '<<lang2(v)<<' '<<1./3<<' '<<0<<"   "<<1/(1+eta)<<' '<<0<<' '<<v<<' '<<(v==0?1./3:eta/(v*(1+eta)))<<' '<<eta<<std::endl;
		if(v<=4) v += dv0; else v *= dv1;
	}
	//         M       eta     eta2    M2         Meta      rho      p    lmbda  Upsilon zeta
	std::cout<<0<<' '<<1<<' '<<1<<' '<<1./3<<' '<<0<<"   "<<.5<<' '<<0<<" inf "<<0<<' '<<1<<std::endl;
	std::cout<<std::endl;

	for(double M4=dM; M4<=.95; M4+=dM){ // .01
		double sp = calc_p(M4); Vec<3> sP(0., 0., sp);
		for(double v=0; v<=20;){  // при v==0 (в приближении среднего поля) известны асимптотики
			double Z=0,  M=0, eta=0, eta2=0, Meta=0, M2=0;
#pragma omp parallel for reduction(+:Z,M,eta,eta2,Meta,M2)
			for(int i=0; i<sz; i++){
				const Vec<3> &m_i = sph_vert(i, rank); double dm_i = sph_vert_area(i, rank);
				for(int j=0; j<sz; j++){ 
					const Vec<3> &m_j = sph_vert(j, rank); double dm_j = sph_vert_area(j, rank);
					double df = exp(sp*(m_i[2]+m_j[2]) + v*m_i*m_j)*dm_i*dm_j;
					Z += df;
					M += m_i[2]*df;
					double e = m_i*m_j;
					eta += e*df;
					eta2 += e*e*df;
					Meta += e*m_i[2]*df;
					M2 += m_i[2]*m_i[2]*df;
				}				
			}
			M /= Z; double p = calc_p(M);
			//         M       eta     eta2    M2                Meta       rho          p        lmbda  Upsilon
			std::cout<<M<<' '<<eta/Z<<' '<<eta2/Z<<' '<<M2/Z<<' '<<Meta/Z<<"   "<<sp/p<<' '<<p<<' '<<v<<' '<<(v==0 ? M/p : (1-sp/p)/v)<<' '<<(eta/Z-M*M)/(1-M*M)<<std::endl;
			if(v<=4) v += dv0; else v *= dv1;
		}
		double M = lang(2*sp);
		//         M       eta     eta2    M2                Meta       rho          p        lmbda  Upsilon zeta
		std::cout<<M<<' '<<1<<' '<<1<<' '<<1-2*M/(2*sp)<<' '<<M<<"   "<<.5<<' '<<(2*sp)<<" inf "<<0<<' '<<1<<std::endl;
		printf("\n");
	}

	// rho неверно!
	for(double v=0; v<=20;){  // при v==0 (в приближении среднего поля) известны асимптотики
		//         M       eta     eta2    M2      Meta       rho  p    lmbda   Upsilon zeta
		std::cout<<1<<' '<<1<<' '<<1<<' '<<1<<' '<<1<<"   "<<.5<<" inf "<<v<<' '<<0<<' '<<1<<std::endl;
		if(v<=4) v += dv0; else v *= dv1;
	}
	//         M       eta     eta2    M2      Meta       rho  p    lmbda   Upsilon
	std::cout<<1<<' '<<1<<' '<<1<<' '<<1<<' '<<1<<"   "<<.5<<" inf inf "<<' '<<0<<' '<<1<<std::endl;
	
	printf("#runtime %g min\n", (omp_get_wtime()-t0)/60);	
	return 0;
}
