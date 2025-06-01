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
	template <typename T> T calc_Mu(T M){
		T M2 = M*M, tmp = 0.775383f + (0.109185f + (0.289114f + (-0.871214f + (1.85968f + (-1.71306f + 0.550912f*M2 )*M2 )*M2 )*M2 )*M2 )*M2;
		return (1-M2*tmp*tmp)/3.;
	}
	template <typename T> T calc_Nu(T M){
		T M2 = M*M;
		return 0.6f + (0.207788f + (0.0567754f + (0.754031f + (-2.99459f + (7.017465f + (-8.202352f + (4.461292f - 0.90036f*M2 )*M2 )*M2 )*M2 )*M2 )*M2 )*M2 )*M2;
	}

	template <typename T> T calc_diff_Mu_div_M(T M){
		T M2 = M*M;
		return -(2.f/3.f)*( 0.775383f + (0.109185f + (0.289114f + (-0.871214f + (1.85968f + (-1.71306f + 0.550912*M2 )*M2 )*M2 )*M2 )*M2 )*M2
						  )*( 0.775383f + (0.327555f + (1.44557f + (-6.098498f  + (16.73712f + ( -18.84366f + 7.161856f*M2 )*M2 )*M2 )*M2 )*M2 )*M2 );
	}
	template <typename T> T calc_diff_Mu(T M){ return M*calc_diff_Mu_div_M<T>(M); }

	template <typename T> T Upsilon(T M, T zeta){
		T M2 = M*M, zeta2 = zeta*zeta, Mz3 = M+zeta; Mz3 *= Mz3*Mz3;
		T delta = 0.1223f*Mz3*Mz3*(1.f-1.6f*M-2.49f*zeta+0.71f*M2+1.71f*zeta2+2.47f*M*zeta-0.12f*M2*M-0.21f*zeta2*zeta-0.34f*M2*zeta-1.17f*M*zeta2);
		return 3*calc_Mu<T>(M)*calc_Mu<T>(zeta)/(1+zeta)*(1+delta);
	}
	template <typename T> T Upsilon(T M, T zeta, T zeta_c){
		T M2 = M*M;
		// return Upsilon<T>(M, zeta)*(1+(0.0256f+(0.0376f-0.1076f*M2)*M2)*zeta*(1-zeta)/(zeta_c*(1-zeta_c)));
		return Upsilon<T>(M, zeta)*(1+(0.0228+(0.0851+(-0.2826+0.1705*M2)*M2)*M2)*zeta*(1-zeta)/(zeta_c*(1-zeta_c)));
	}
	template <typename T> T Upsilon(T M, T zeta, T zeta_c, T eta, T eta_c){
		T U = Upsilon<T>(M, zeta, zeta_c);
		/* if(zeta>zeta_c){
			T Mc2 = (eta-eta_c)/(1-eta_c);
			if(Mc2>0){
				T Mc = sqrtf(Mc2), m = 1.f - M/Mc, e1e4 = eta*(1-eta); e1e4 *= e1e4; e1e4 *= e1e4;
				U += 0.714f*e1e4*m*(1-3.2686f*m+8.1953f*eta+15.7007f*m*eta);
			}
		}*/
		return U;
	}
	template <typename T> T Q_MFA(T M, int n_b){ T mu_p = calc_Mu<T>(M); return mu_p*(2 + (n_b-1)*M*M - 3*mu_p); }
	template <typename T> T Qiji(T M, T zeta, T mu_p, T mu_z){
		T M2 = M*M, zeta2 = zeta*zeta;
		return .5f*mu_p*mu_z*(6*(2-3*mu_p) -10.16f*zeta*M2*(1-1.217f*M-1.196f*zeta+1.929f*M*zeta+1.329f*M2+0.836f*zeta2-1.077f*M2*zeta-0.761f*M*zeta2));
	}
	//--------------------------------------------------------------------------
	//  single axis case
	//--------------------------------------------------------------------------
	template <typename T> struct SingleAxisCoeffLLB {
		T M;           ///< input parametr <m>
		T mu_p, nu_p;  ///< middle output parameters
		T Theta, Xi;   ///< LLB output parameters (1st equation)
		
		void calc(){
			mu_p = calc_Mu<T>(M); nu_p = calc_Nu<T>(M); 
			Theta = 2*M*mu_p*nu_p; Xi = 2*mu_p;
		}
	};
	//--------------------------------------------------------------------------
	template <typename T> struct BaseCoeffCMD {
		int n_b; T eta_c;  ///< input parametrs of material 
		T eta;             ///< input parametrs
		T q;   ///< pankrat coeff
		
		T zeta, zeta_c;    ///< middle output parameters
		T U, Q;            ///< CMD output parameters
	protected:
		void calc(T M, T M2, T mu_p){		
			zeta = (eta-M2)/(1-M2); zeta_c = eta_c -.135f*sqrtf(fabsf(M))*M +.123f*M2*M;

			U = Upsilon<T>(M, zeta, zeta_c, eta, eta_c);
			if(zeta<=zeta_c){			
				T Q0 = Q_MFA<T>(M, n_b), ec = M2 + (1-M2)*zeta_c,  Uc = Upsilon<T>(M, zeta_c, zeta_c, ec, eta_c), Q1 = n_b*ec*Uc;
				
				T dzeta = 1e-4f, e1 = M2 + (1-M2)*(zeta_c+dzeta), z1 = (e1-M2)/(1-M2);
				T dQdz = (n_b*e1*Upsilon(M, z1, zeta_c, e1, eta_c) - Q1)/dzeta; // тут все таки нужна аналитическая производная
					
				T B = (Q0+dQdz*zeta_c-Q1)/(zeta_c*zeta_c), A = dQdz - 2*B*zeta_c;   Q = Q0 + (A + B*zeta)*zeta;
				// T B = (Q0+dQdz*zeta_c-Q1)/(2*zeta_c*zeta_c*zeta_c), A = dQdz - 3*B*zeta_c*zeta_c;	Q = Q0 + (A + B*zeta*zeta)*zeta;
			} else {
				// T Mc2 = (eta-eta_c)/(1-eta_c);
				Q = n_b*eta*U; // - .1*(Mc2-M2)*(Mc2-1.25*M2)*(1-eta); // *(eta-eta_c)/(1-eta_c); // *(1 - 0.2*(1-zeta)*(zeta-zeta_c)*(1-sqrt(fabs(M))));
				/*
				T Mc2 = (eta-eta_c)/(1-eta_c);
				if(Mc2>0){
					T Mc = sqrtf(Mc2), m = 1.f - M/Mc, eta2 = eta*eta, m2 = m*m;
					// Q += -0.2186f*eta2*(1-eta)*m*(1+1.3699f*m+1.1961f*eta-10.3134f*m*eta+4.6053f*m2+3.9136f*eta2) + .01f*m;
					Q += .01f*m;
				}
				*/
				
				
				/*
				T Mc = sqrtf((eta-eta_c)/(1-eta_c));  // eta_c = (eta-Mc2)/(1-Mc2)  ==>  eta_c - Mc2*eta_c = eta - Mc2  ==>  Mc2 = (eta-eta_c)/(1-eta_c)
				T a = 0, b = .999;  // if(a<0){ a = 0; } if(b>1){ b = 1; }
				for(int i=0; i<10; i++){
					T zc = eta_c -.135f*sqrtf(fabsf(Mc))*Mc +.123f*Mc*Mc*Mc;
					if(zc<(eta-Mc*Mc)/(1-Mc*Mc)) a = Mc; else b = Mc;
					Mc = (a+b)/2;
				}
				T Mc2 = Mc*Mc, zc = (eta-Mc2)/(1-Mc2), Uc = Upsilon<T>(Mc, zc, zc);
				Q = n_b*eta*Uc - .02f*(1-M2/Mc2); // *(1 - 0.4*(1-zeta)*(zeta-zeta_c)*(Mc-M)); // n_b*eta*(U*.25+Uc*.75);
				*/
				// U = Uc;
			}			
			/*
			T dUM = (0.0256f+(0.0376f-0.1076f*M2)*M2), dUMz = dUM/(zeta_c*(1-zeta_c)); // *zeta_c*(1-zeta_c));
			// T dUM = 0, dUMz = dUM/(zeta_c*(1-zeta_c)); // *zeta_c*(1-zeta_c));
			// T wM = 0, dUMz2 = dUMz/(zeta_c*(1-zeta_c))*wM;  dUMz *= 1-wM; 
			// U = 3*mu_p*calc_Mu<T>(zeta)/(1+zeta)*(1+(dUMz + dUMz2*zeta*(1-zeta))*zeta*(1-zeta)); //   + q*zeta*(zeta-zeta_c)*(1-zeta); // *zeta*(1-zeta));			
			U = 3*mu_p*calc_Mu<T>(zeta)/(1+zeta)*(1+dUMz*zeta*(1-zeta)); //   + q*zeta*(zeta-zeta_c)*(1-zeta); // *zeta*(1-zeta));			

			// V0
			if(zeta<=zeta_c){			
				T ec = M2 + (1-M2)*zeta_c, mu_zc = calc_Mu<T>(zeta_c), _1zc = 1/(1+zeta_c); 
				T Q0 = Q_MFA<T>(M, n_b), Uc = 3*mu_p*mu_zc*_1zc*(1+dUM), Q1 = n_b*ec*Uc;
				
				// T dzeta = 1e-4f, e1 = M2 + (1-M2)*(zeta_c+dzeta), z1 = (e1-M2)/(1-M2);
				// T dQdz = (n_b*e1*(3*mu_p*calc_Mu<T>(z1)/(1+z1)*(1+(dUMz+dUMz2*z1*(1-z1))*z1*(1-z1))) - Q1)/dzeta; // тут все таки нужна аналитическая производная
				T dQdz = n_b*((1-M2)*Uc + ec*3*mu_p*_1zc*( (calc_diff_Mu(zeta_c)-mu_zc*_1zc)*(1+dUM) + mu_zc*dUMz*(1-2*zeta_c) ));
					
				T B = (Q0+dQdz*zeta_c-Q1)/(zeta_c*zeta_c), A = dQdz - 2*B*zeta_c;   Q = Q0 + (A + B*zeta)*zeta;
				// T B = (Q0+dQdz*zeta_c-Q1)/(2*zeta_c*zeta_c*zeta_c), A = dQdz - 3*B*zeta_c*zeta_c;	Q = Q0 + (A + B*zeta*zeta)*zeta;
			} else {  // Q = n_b*eta*U; // *(1 - 0.2*(1-zeta)*(zeta-zeta_c)*(1-sqrt(fabs(M))));

				
				T Mc = (eta-eta_c)/(1-eta_c);  // eta_c = (eta-Mc2)/(1-Mc2)  ==>  eta_c - Mc2*eta_c = eta - Mc2  ==>  Mc2 = (eta-eta_c)/(1-eta_c)
				T a = 0, b = .999;  // if(a<0){ a = 0; } if(b>1){ b = 1; }
				for(int i=0; i<10; i++){
					T zc = eta_c -.135f*sqrtf(fabsf(Mc))*Mc +.123f*Mc*Mc*Mc;
					if(zc<(eta-Mc*Mc)/(1-Mc*Mc)) a = Mc; else b = Mc;
					Mc = (a+b)/2;
				}
				T Mc2 = Mc*Mc, zc = (eta-Mc2)/(1-Mc2), mu_zc = calc_Mu<T>(zc), _1zc = 1/(1+zc); 
				T Uc = 3*calc_Mu<T>(Mc)*mu_zc*_1zc*(1+(0.0256f+(0.0376f-0.1076f*Mc2)*Mc2));
				Q = n_b*eta*Uc; // *(1 - 0.4*(1-zeta)*(zeta-zeta_c)*(Mc-M)); // n_b*eta*(U*.25+Uc*.75);
				
				U = Uc;
			}
			*/
			// Q *= 1 - .1*sqrt(fabs(M))*zeta*(zeta-zeta_c)*(1-zeta)-.2*M*M*(zeta>zeta_c)*(zeta-zeta_c)*(zeta-zeta_c)*(1-zeta);
			// if(zeta>zeta_c) U *= 1-(zeta-zeta_c)*(1-zeta)*(1-M);
			// Q *= 1 - sqrt(fabs(M))*zeta*(zeta-zeta_c)*(1-zeta);

			// if(zeta<zeta_c) Q *= 1 - .25*sqrt(fabs(M))*zeta*(zeta-zeta_c)*(1-zeta);
			// else Q *= 1 + q*(zeta-zeta_c)*(zeta-zeta_c)*(1-eta*eta*eta/2); //sqrt(fabs(M)));

			/* 
			else{
				// T sc = eta>eta_c*1.45f ? std::max(-2.,-.05/(eta-eta_c*1.45)): -2.;  // BCC
				T sc = eta>eta_c*1.5f ? std::max(-2.,-.05/(eta-eta_c*1.5)): -2.;  // BCC
				// T sc = eta>eta_c*1.2 ? std::max(-2.,-.025/(eta-eta_c*1.2)): -2.;
				Q *= 1 +  sc*(zeta-zeta_c)*(zeta-zeta_c)*(1-eta*eta*eta/2); //sqrt(fabs(M)));
			} */
			// Q *= 1 - .25*pow(fabs(M), .25)*zeta*(zeta-zeta_c)*(1-zeta);
			// if(zeta<zeta_c) Q *= 1 - pow(M, .25)*zeta*(zeta-zeta_c);


			/*  бол-мен правильные зависимости при H,K=0 и местами правильная восприимчивость но есть артефакты
			// V1
			// U = 3*mu_p*calc_Mu<T>(zeta)/(1+zeta); // *zeta*(1-zeta));			
			T mu_z = calc_Mu<T>(zeta),  mu_zc = calc_Mu<T>(zeta_c), ec = M2 + (1-M2)*zeta_c;
			T Q0 = Qiji(M, zeta, mu_p, mu_z), Q0c = Qiji(M, zeta_c, mu_p, mu_zc);
			T Qan = mu_p*M2*(1-zeta)*(1-zeta), Qanc = mu_p*M2*(1-zeta_c)*(1-zeta_c);
			T Qtr = calc_Mu<T>(sqrt(eta))*eta, Qtrc = calc_Mu<T>(sqrt(ec))*ec;
			T Qc = n_b*ec*3*mu_p*mu_zc/(1+zeta_c); // *(1+dUM);
			int nb1 = n_b-1;
			if(n_b==12){ Qc -= 4*Qtrc; nb1 = 7; }
			// (n_b-1)*((1-w)*Qan + w*Qtr) = Q-Q0 --> w = ((Q-Q0)/(n_b-1)-Qan)/(Qtr-Qan)    
			T w = ((Qc-Q0c)/nb1-Qanc)/(Qtrc-Qanc) * zeta/zeta_c;
			Q = Q0 + nb1*((1-w)*Qan + w*Qtr);
			if(n_b==12) Q += 4*Qtr;
			*/
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
			mu_p = calc_Mu<T>(Mabs); nu_p = calc_Nu<T>(Mabs);
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
