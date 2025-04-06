#ifndef CMDLIB_CMD_COEFF_HPP
#define CMDLIB_CMD_COEFF_HPP
/**
 * Copyright (C) 2025 Antov V. Ivanov  <aiv.racs@gmail.com>
 * Licensed under the Apache License, Version 2.0
 *
 * Header-only aiwlib-free implementation of CMD coefficients.
 **/

#include <math.h>

namespace CMD{
	//--------------------------------------------------------------------------
	//  vector operators
	//--------------------------------------------------------------------------
	template <typename T> T sprod(const T* a, const T* b){ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
	template <typename T> void vprod(const T* a, const T* b, T* c){
		c[0] = a[1]*b[2]-a[2]*b[1];
		c[1] = a[2]*b[0]-a[0]*b[2];
		c[2] = a[0]*b[1]-a[1]*b[0];
	}
	//--------------------------------------------------------------------------
	//  special functions approximation
	//--------------------------------------------------------------------------
	template <typename T> T calc_mu_p(T M){
		T M2 = M*M, tmp = 0.775383f + (0.109185f + (0.289114f + (-0.871214f + (1.85968f + (-1.71306f + 0.550912f*M2 )*M2 )*M2 )*M2 )*M2 )*M2;
		return (1-M2*tmp*tmp)/3.;
	}
	template <typename T> T calc_nu_p(T M){
		T M2 = M*M;
		return 0.6f + (0.207788f + (0.0567754f + (0.754031f + (-2.99459f + (7.017465f + (-8.202352f + (4.461292f - 0.90036f*M2 )*M2 )*M2 )*M2 )*M2 )*M2 )*M2 )*M2;
	}

	template <typename T> T calc_diff_mu_p_div_M(T M){
		T M2 = M*M;
		return -(2.f/3.f)*( 0.775383f + (0.109185f + (0.289114f + (-0.871214f + (1.85968f + (-1.71306f + 0.550912*M2 )*M2 )*M2 )*M2 )*M2 )*M2
						  )*( 0.775383f + (0.327555f + (1.44557f + (-6.098498f  + (16.73712f + ( -18.84366f + 7.161856f*M2 )*M2 )*M2 )*M2 )*M2 )*M2 );
	}
	template <typename T>  T calc_diff_mu_p(T M){ return M*calc_diff_mu_p_div_M<T>(M); }

	template <typename T> T Upsilon(T M, T zeta){ return 3*calc_mu_p<T>(M)*calc_mu_p<T>(zeta)/(1+zeta); }
	template <typename T> T Upsilon(T M, T zeta, T zeta_c){
		T M2 = M*M;
		return Upsilon<T>(M, zeta)*(1+(0.0256f+(0.0376f-0.1076f*M2)*M2)*zeta*(1-zeta)/(zeta_c*(1-zeta_c)));
	}
	template <typename T> T Q_MFA(T M, int n_b){ T mu_p = calc_mu_p<T>(M); return mu_p*(2 + (n_b-1)*M*M - 3*mu_p); }
	//--------------------------------------------------------------------------
	//  single axis case
	//--------------------------------------------------------------------------
	template <typename T> struct SingleAxisCoeffLLB {
		T M;           ///< input parametr <m>
		T mu_p, nu_p;  ///< middle output parameters
		T Theta, Xi;   ///< LLB output parameters (1st equation)
		
		void calc(){
			mu_p = calc_mu_p<T>(M); nu_p = calc_nu_p<T>(M); 
			Theta = 2*M*mu_p*nu_p; Xi = 2*mu_p;
		}
	};
	//--------------------------------------------------------------------------
	template <typename T> struct BaseCoeffCMD {
		int n_b; T eta_c;  ///< input parametrs of material 
		T eta;             ///< input parametrs

		T zeta, zeta_c;    ///< middle output parameters
		T U, Q;            ///< CMD output parameters
	protected:
		void calc(T M, T M2, T mu_p){		
			zeta = (eta-M2)/(1-M2); zeta_c = eta_c -.135f*sqrtf(M)*M +.123f*M2*M;
			T dUM = 0.0256f+(0.0376f-0.1076f*M2)*M2, dUMz = dUM/(zeta_c*(1-zeta_c));
			U = 3*mu_p*calc_mu_p<T>(zeta)/(1+zeta)*(1+dUMz*zeta*(1-zeta));			

			if(zeta<=zeta_c){			
				T ec = M2 + (1-M2)*zeta_c, mu_zc = calc_mu_p<T>(zeta_c), _1zc = 1/(1+zeta_c); 
				T Q0 = Q_MFA<T>(M, n_b), Uc = 3*mu_p*mu_zc*_1zc*(1+dUM), Q1 = n_b*ec*Uc;
				// T dzeta = 1e-4f, e1 = M2 + (1-M2)*(zeta_c+dzeta), z1 = (e1-M2)/(1-M2);
				// T dQdz = (n_b*e1*3*mu_p*calc_mu_p<T>(z1)/(1+z1)*(1+dUMz*z1*(1-z1)) - Q1)/dzeta; // тут все таки нужна аналитическая производная
				T dQdz = n_b*((1-M2)*Uc + ec*3*mu_p*_1zc*( (calc_diff_mu_p(zeta_c)-mu_zc*_1zc)*(1+dUM) + mu_zc*dUMz*(1-2*zeta_c) ));
					
				T B = (Q0+dQdz*zeta_c-Q1)/(zeta_c*zeta_c), A = dQdz - 2*B*zeta_c;
				Q = Q0 + (A + B*zeta)*zeta;
			} else Q = n_b*eta*U;
		}
	};	
	//--------------------------------------------------------------------------
	template <typename T> struct SingleAxisCoeffCMD: public SingleAxisCoeffLLB<T>, public BaseCoeffCMD<T> {
		using SingleAxisCoeffLLB<T>::M;
		using SingleAxisCoeffLLB<T>::mu_p;
		using BaseCoeffCMD<T>::eta;
		
		T Psi;       ///< CMD output parameters
		
		void calc(){  // вынести это в отдельную шаблонную функцию что бы не было дублирования?
			SingleAxisCoeffLLB<T>::calc();
			T M2 = M*M;  BaseCoeffCMD<T>::calc(M, M2, mu_p);
			Psi = 0.9226f*(1-eta)*M2*M;
		}
	};	
	//--------------------------------------------------------------------------
	//   3D case
	//--------------------------------------------------------------------------
	template <typename T>  struct CoeffLLB {
		T M[3], H[3], nK[3];   ///< input parametrs

		T M2, MnK, Mabs, mu_p, nu_p;  ///< middle output parameters
		T PHI[3], THETA[3], XiH[3];   ///< LLB output parameters (1st equation)

		void calc(){
			M2 = 0; for(int i=0; i<3; i++){ M2 += M[i]*M[i]; }  Mabs = sqrt(M2);
			mu_p = calc_mu_p<T>(Mabs); nu_p = calc_nu_p<T>(Mabs);
			MnK = sprod(M, nK); T MMnK[3]; vprod(M, nK, PHI); vprod(M, PHI, MMnK);
			T t = MnK*nu_p; for(int i=0; i<3; i++) PHI[i] *= t;
			if(M2>1e-6f){
				t = MnK/M2; T mnp = mu_p*nu_p, t1 = (5*t*MnK-1)*mnp, t2 = 2*mnp*MnK;
				for(int i=0; i<3; i++) THETA[i] = t1*M[i] - t2*nK[i] - t*MMnK[i];
			} else for(int i=0; i<3; i++) THETA[i] = 0;
			for(int i=0; i<3; i++){
				XiH[i] = H[i]*(1-mu_p);
				for(int j=0; j<3; j++) XiH[i] -= H[j]*M[i]*M[j]*nu_p;
			}
		}
	};
	//--------------------------------------------------------------------------
	template <typename T>  struct CoeffCMD: public CoeffLLB<T>,  public BaseCoeffCMD<T> {
		using CoeffLLB<T>::Mabs;
		using CoeffLLB<T>::M2;
		using CoeffLLB<T>::MnK;
		using CoeffLLB<T>::mu_p;
		using BaseCoeffCMD<T>::eta;
		
		T Psi;     ///< CMD output parameters

		void calc(){
			CoeffLLB<T>::calc();			
			Psi = (1.3836f*MnK*MnK - .46134f*M2)*(1-eta)*Mabs;
			BaseCoeffCMD<T>::calc(Mabs, M2, mu_p);
		}
	};
	//--------------------------------------------------------------------------
}

#endif // CMDLIB_CMD_COEFF_HPP
