#include "../include/llb_coeff.hpp"
using namespace aiw;
//------------------------------------------------------------------------------
// const float MAX_P = logf(HUGE_VALF)-1;
double calc_p(double M){
	if(fabs(M)<1e-6) return M*3;
	double p = 1.f;
	for(int i=0; i<10; i++){
		double Lp = calc_L(p);
		if(p>400 || fabs(Lp-M)<1e-6) break;
		p = p-(Lp-M)/calc_dLdp(p);
		if(fabs(p)<1e-6) return p;
	}
	return p;
}
//------------------------------------------------------------------------------
const float _3 = 1./3;
float calc_mu_p(float M2){
	float tmp = 0.775383f + (0.109185f + (0.289114f + (-0.871214f + (1.85968f + (-1.71306f + 0.550912f*M2 )*M2 )*M2 )*M2 )*M2 )*M2;
	return (1-M2*tmp*tmp)*_3;
}
float calc_nu_p(float M2){
	return 0.6f + (0.207788f + (0.0567754f + (0.754031f + (-2.99459f + (7.017465f + (-8.202352f + (4.461292f - 0.90036f*M2 )*M2 )*M2 )*M2 )*M2 )*M2 )*M2 )*M2;
}
//------------------------------------------------------------------------------
int LLBCoeff::n_b; 
aiw::Vecf<3> LLBCoeff::nK; 
//------------------------------------------------------------------------------
void LLBCoeff::init(float M_, bool exact){
	M = fabsf(M_); float M2 = M*M;
	if(exact){
		double p = calc_p(M);
		mu_p = M>1e-6? M/p: calc_mu_p(M2);
		nu_p = M>1e-6? (1-3*mu_p)/(M*M): calc_nu_p(M2);
		S = log(4*M_PI*sinh(p)) - p*M;
	} else {
		mu_p = calc_mu_p(M2);
		nu_p = calc_nu_p(M2);
		S = log(4*M_PI*sinh(M/mu_p)) - M2/mu_p;  // аппроксимация?
	}
	Phi = Vecf<1>();
	Theta = vecf(0, 0, 2*M*mu_p*nu_p);
	XI.fill(); XI(2, 2) = 2*mu_p;

	Q = -mu_p*(3*mu_p - (n_b-1)*M2 - 2);
	// float S;
}
//------------------------------------------------------------------------------
void LLBCoeff::init(const aiw::Vecf<3> &M_, bool exact){
	float M2 = M_*M_, MnK = M_*nK; M = sqrt(M2); 
	if(exact){
		double p = calc_p(M);
		mu_p = M>1e-6? M/p: calc_mu_p(M2);
		nu_p = M>1e-6? (1-3*mu_p)/(M*M): calc_nu_p(M2);
		S = log(4*M_PI*sinh(p)) - p*M;
	} else {
		mu_p = calc_mu_p(M2);
		nu_p = calc_nu_p(M2);
		S = log(4*M_PI*sinh(M/mu_p)) - M2/mu_p;  // аппроксимация?
	}
	Phi = MnK*nu_p*(M_%nK);
	if(M>1e-6) Theta = mu_p*nu_p*(5*MnK*MnK/M2*M_ - M_ - 2*MnK*nK) - M_%(M_%nK)*(MnK/M2);
	else Theta = vecf(0.f); // ???
	
	for(int i=0; i<3; i++){
		XI(i, i) = 1-mu_p - M_[i]*M_[i]*nu_p;
		for(int j=0; j<i; j++) XI(i, j) = XI(j, i) = -M_[i]*M_[j]*nu_p;
	}
	
	Q = -mu_p*(3*mu_p - (n_b-1)*M2 - 2);
}
//------------------------------------------------------------------------------

