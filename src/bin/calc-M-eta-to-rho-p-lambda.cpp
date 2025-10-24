#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string>
#include <iostream>
#include <aiwlib/sphere>

#include "include/Z2.hpp"

inline double Z(double p){ return 4*M_PI*(fabs(p)>1e-4? sinh(p)/p: 1+p*p/6); }

double calc_mu_p(double M){
	double M2 = M*M, tmp = 0.775383f + (0.109185f + (0.289114f + (-0.871214f + (1.85968f + (-1.71306f + 0.550912f*M2 )*M2 )*M2 )*M2 )*M2 )*M2;
	return (1-M2*tmp*tmp)/3.;
}

int main(int argc, const char **argv){
	if(argc!=3 && (argc!=4 || std::string(argv[3])!="angle")){
		std::cerr<<"usage: ./calc-M-eta-to-rho-p-lambda M eta [angle] ==> rho p lambda Z2 S1 S2 M eta eta2  [rho p lambda Z2 S1 S2 M eta eta2 (for angle)]\n";
		exit(1);
	}
	double M = atof(argv[1]), eta = atof(argv[2]), p = invL(M), zeta = (eta-M*M)/(1-M*M);
	Z2_symmetrical z2;  z2.m = M; z2.eta = eta;
	double rho = 1-zeta/2-0.514*(1-zeta)*zeta*( 1 - 0.232*M + 1.681*M*M - 1.466*M*M*M - 0.123*M*zeta
												- 1.041*zeta + 0.858*zeta*zeta - 0.356*zeta*zeta*zeta);
	
	z2.calc_from_moments(p*rho, (1-rho)/(3*calc_mu_p(M)*calc_mu_p(zeta)/(1+zeta))); double z = z2.Z2_norm/exp(-2*z2.h-z2.l), S2 = log(z)/2 - z2.h*M -z2.l*eta/2;
	printf("{'rho':%f, 'p':%f, 'l':%f, 'Z2':%f, 'S1':%f, 'S2':%f, 'eta2':%f", (p>1e-3? z2.h/p: 1), p, (eta>M*M?z2.l:0), z, log(Z(p))-p*M, S2, z2.eta2);

	if(argc==4 && std::string(argv[3])=="angle"){
		// f3_angle = exp( rho p m_i + rho p m_k + (2rho-1) p m_j + lambda m_i m_j + lambda m_j m_k) / Z2^2 * Z1
		// f2_ik = exp( rho p m_i + rho p m_k) / Z2^2 * Z1 * Z(| (2rho-1) p + lambda m_i + lambda m_k |)
		int rank = 3;  rho = p>1e-3? z2.h/p: 1;
		aiw::sph_init_table(rank);  int sz = aiw::sph_vertex_num(rank);  M = eta = z = 0;
#pragma omp parallel for reduction(+:z,M,eta)
		for(int i=0; i<sz; i++){
			const aiw::Vec<3> &m_i = aiw::sph_vert(i, rank); double dm_i = aiw::sph_vert_area(i, rank);
			for(int j=0; j<sz; j++){
				const aiw::Vec<3> &m_j = aiw::sph_vert(j, rank); double dm_j = aiw::sph_vert_area(j, rank);
				double df = exp(z2.h*(m_i[2]+m_j[2]))*Z((aiw::vec(0., 0., (2*rho-1)*p)+z2.l*(m_i+m_j)).abs())*dm_i*dm_j;
				z += df;
				M += m_i[2]*df;
				eta += m_i*m_j*df;
			}
		}
		M /= z; eta /= z; p = invL(M), zeta = (eta-M*M)/(1-M*M);
		z2.m = M; z2.eta = eta;
		rho = 1-zeta/2-0.514*(1-zeta)*zeta*( 1 - 0.232*M + 1.681*M*M - 1.466*M*M*M - 0.123*M*zeta
											 - 1.041*zeta + 0.858*zeta*zeta - 0.356*zeta*zeta*zeta);	
		z2.calc_from_moments(p*rho, (1-rho)/(3*calc_mu_p(M)*calc_mu_p(zeta)/(1+zeta))); z = z2.Z2_norm/exp(-2*z2.h-z2.l); S2 = log(z)/2 - z2.h*M -z2.l*eta/2;
		printf(", 'rho_a':%f, 'l_a':%f, 'Z2a':%f, 'S2a':%f, 'eta_a':%f, 'eta2a':%f", (p>1e-3? z2.h/p: 1), (eta>M*M?z2.l:0), z, S2, z2.eta, z2.eta2);
	}
	printf("}\n");
	
	return 0;
}
