#ifndef CMDLIB_LLB_COEFF_HPP
#define CMDLIB_LLB_COEFF_HPP

/**
 * Copyright (C) 2025 Antov V. Ivanov  <aiv.racs@gmail.com>
 * Licensed under the Apache License, Version 2.0
 **/

#include <aiwlib/vec>
#include <aiwlib/matr>
/*
нужно две версии - точная и приблизительная, с одинаковым интерфейсом
отличаться они могут только вычислением p
варианты - функции принимающие f1Exact/f1Approx ?
либо структура с параметризованным 
 */
//------------------------------------------------------------------------------
// -def mu_p=(1-(0.775383*M+0.109185*M**3+0.289114*M**5-0.871214*M**7+1.85968*M**9-1.71306*M**11+0.550912*M**13)**2)/3.
// -def diff_mu_p=-2./3.*(0.775383*M+0.109185*M**3+0.289114*M**5-0.871214*M**7+1.85968*M**9-1.71306*M**11+0.550912*M**13)*(7.161856*M**12-18.84366*M**10+16.73712*M**8-6.098498*M**6+1.44557*M**4+0.327555*M**2+0.775383)
// -def Q_MFA=mu_p*(3*mu_p-(n_b-1)*M**2-2)
//------------------------------------------------------------------------------
inline double calc_Z(double p){ return 4*M_PI*(fabs(p)>1e-4? sinh(p)/p: 1+p*p/6); }
inline double calc_L(double p){ return fabs(p)>1e-4? 1/tanh(p)-1/p: p/3; }
inline double calc_dLdp(double p){ double shp = sinh(p); return fabs(p)>1e-6? 1/(p*p)-1/(shp*shp) : 1/3.-p*p/15; }
double calc_p(double M);
float calc_mu_p(float M2);
float calc_nu_p(float M2);
//------------------------------------------------------------------------------
struct LLBCoeff {
	static int n_b; 
	static aiw::Vecf<3> nK;  ///< input parametr n_K

	float M; 
	float mu_p;
	float nu_p;

	aiw::Vecf<3> Phi;
	aiw::Vecf<3> Theta;
	aiw::Matr<3, 3, float> XI;
	float Q;
	float S;
	
	void init(float M_, bool exact=false);
	void init(const aiw::Vecf<3> &M_, bool exact=false);
};
//------------------------------------------------------------------------------
#endif // CMDLIB_LLB_COEFF_HPP
