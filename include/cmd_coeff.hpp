#ifndef CMDLIB_CMD_COEFF_HPP
#define CMDLIB_CMD_COEFF_HPP
/**
 * Copyright (C) 2025 Antov V. Ivanov  <aiv.racs@gmail.com>
 * Licensed under the Apache License, Version 2.0
 *
 * Header-only aiwlib-free implementation of CMD coefficients.
 **/

#include <math.h>
// #include <aiwlib/debug>

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
		T mu_p, nu_p;  ///< middle output parameters
		T Theta, Xi;   ///< LLB output parameters (1st equation)
		
		void calc(T M){  ///< M --- input parametr <m>
			mu_p = calc_Mu<T>(M); nu_p = calc_Nu<T>(M); 
			Theta = 2*M*mu_p*nu_p; Xi = 2*mu_p;
		}
	};
	//--------------------------------------------------------------------------
	template <typename T> struct ConstCMD {
		int n_b;              ///< number of nearest neighbors
		T eta_c;              ///< input parametrs of material
		T nK[3] = {0, 0, 1};  ///< anisotropy axis
		T u = 0, q = 0;       ///< corrections in the nonequilibrium region
	};
	//--------------------------------------------------------------------------
	template <typename T> struct BaseCoeffCMD {
		T zeta, zeta_c;    ///< middle output parameters
		T U, Q;            ///< CMD output parameters
	protected:
		void calc(T M, T eta, T M2, T mu_p, const ConstCMD<T>& cell){		
			zeta = (eta-M2)/(1-M2); zeta_c = cell.eta_c -.135f*sqrtf(fabsf(M))*M +.123f*M2*M;

			U = Upsilon<T>(M, zeta, zeta_c, eta, cell.eta_c);
			if(zeta<=zeta_c){			
				T Q0 = Q_MFA<T>(M, cell.n_b), ec = M2 + (1-M2)*zeta_c,  Uc = Upsilon<T>(M, zeta_c, zeta_c, ec, cell.eta_c), Q1 = cell.n_b*ec*Uc;
				
				T dzeta = 1e-4f, e1 = M2 + (1-M2)*(zeta_c+dzeta), z1 = (e1-M2)/(1-M2);
				T dQdz = (cell.n_b*e1*Upsilon(M, z1, zeta_c, e1, cell.eta_c) - Q1)/dzeta; // тут все таки нужна аналитическая производная ?
				// T dQdz = cell.n_b*((1-M2)*Uc + ec*3*mu_p*_1zc*( (calc_diff_Mu(zeta_c)-mu_zc*_1zc)*(1+dUM) + mu_zc*dUMz*(1-2*zeta_c) ));
					
				T B = (Q0+dQdz*zeta_c-Q1)/(zeta_c*zeta_c), A = dQdz - 2*B*zeta_c;   Q = Q0 + (A + B*zeta)*zeta;
			} else {
				Q = cell.n_b*eta*U;
				T dM = 1-M, dz = zeta-zeta_c, dz2 = dz*dz;
				U *= 1 + cell.u*dM*dM*dM*dz2;
				Q *= 1 + cell.q*dM*dz2*dz;
			}

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
		using SingleAxisCoeffLLB<T>::mu_p;
		
		T Psi;       ///< CMD output parameters
		
		void calc(T M, T eta, const ConstCMD<T>& cell){  // вынести это в отдельную шаблонную функцию что бы не было дублирования?
			SingleAxisCoeffLLB<T>::calc(M);
			T M2 = M*M;  BaseCoeffCMD<T>::calc(M, eta, M2, mu_p, cell);
			Psi = 0.9226f*(1-eta)*M2*M;
		}
	};	
	//--------------------------------------------------------------------------
	//   3D case
	//--------------------------------------------------------------------------
	template <typename T>  struct CoeffLLB {
		T M2, MnK, Mabs, mu_p, nu_p;  ///< middle output parameters
		T PHI[3], THETA[3], XiH[3];   ///< LLB output parameters (1st equation)

		void calc(const T M[3], const T H[3], const T nK[3]){
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
		
		T Psi;     ///< CMD output parameters

		void calc(const T M[3], T eta, const T H[3], const ConstCMD<T> &cell){
			CoeffLLB<T>::calc(M, H, cell.nK);			
			Psi = (1.3836f*MnK*MnK - .46134f*M2)*(1-eta)*Mabs;
			BaseCoeffCMD<T>::calc(Mabs, eta, M2, mu_p, cell);
		}
	};
	//--------------------------------------------------------------------------
	//   frontend (right sides of CMD equation)
	//--------------------------------------------------------------------------
	template <typename Type> struct Right: public ConstCMD<Type> {
		using ConstCMD<Type>::n_b;
		using ConstCMD<Type>::nK;
		using ConstCMD<Type>::eta_c;
		using ConstCMD<Type>::u;
		using ConstCMD<Type>::q;
		
		Type eG;
		Type alpha;
		Type T;
		Type K;
		Type J = 1;
		Type CdMdt; 
		// gamma = 1 ???
		
		//----------------------------------------------------------------------
		// void init_SC(){  n_b = 6;  eta_c = .33;  eG = .725;  u = 1.74; q = -0.36; CdMdt = -2.29; }
		// void init_BCC(){ n_b = 8;  eta_c = .27;  eG = .775;  u = 2.2;  q =  0.81; CdMdt = -1.72; }
		// void init_FCC(){ n_b = 12; eta_c = .25;  eG = .7975; u = 1.42; q =  0.79; CdMdt = -1.22; }
		void init_SC(){  n_b = 6;  eta_c = .33;  eG = .725;  u = 1.74097; q = -0.360312; CdMdt = -2.29127/2; }
		void init_BCC(){ n_b = 8;  eta_c = .27;  eG = .775;  u = 2.19811; q =  0.813378; CdMdt = -1.71897/2; }
		void init_FCC(){ n_b = 12; eta_c = .25;  eG = .7975; u = 1.41678; q =  0.794254; CdMdt = -1.21721/2; }
		//----------------------------------------------------------------------
		void single_axis_LLB(Type H, Type M, Type &dM) const {
			SingleAxisCoeffLLB<Type> coeff; coeff.calc(M);
			dM = 2*alpha*(K*coeff.Theta + .5f*coeff.Xi*(H + n_b*eG*M) - T*M);
		}
		void single_axis_CMD(Type H, Type M, Type eta, Type &dM, Type &deta) const {
			SingleAxisCoeffCMD<Type> coeff;  coeff.calc(M, eta, (const ConstCMD<Type>&)*this);
			// WMSG(coeff.mu_p, coeff.nu_p, coeff.Theta, coeff.Xi, coeff.Psi, coeff.U, coeff.Q);
			dM   = 2*alpha*(K*coeff.Theta + .5*coeff.Xi*H + (n_b*J*coeff.U - T)*M);   //dMdt*U
			deta = 4*alpha*(coeff.U*H*M + K*coeff.Psi + J*coeff.Q - T*eta);
			if(coeff.zeta>coeff.zeta_c){
				Type C = CdMdt*(coeff.zeta-coeff.zeta_c);
				deta += 4*alpha*C*dM*H;
				dM += 2*alpha*n_b*J*C*dM;
			}
		}
		//----------------------------------------------------------------------
		void LLB(const Type H[3], const Type M[3], Type dM[3]) const {
			CoeffLLB<Type> coeff; coeff.calc(M, H, nK);  Type MxH[3], c = eG*n_b*J*coeff.mu_p; vprod(M, H, MxH);
			for(int i=0; i<3; i++) dM[i] = -MxH[i] + 2*K*(coeff.PHI[i] + alpha*coeff.THETA[i]) + alpha*(coeff.XiH[i] + 2*(c - T)*M[i]); 
		}
		void CMD(const Type H[3], const Type M[3], Type eta, Type dM[3], Type &deta) const {
			CoeffCMD<Type> coeff; coeff.calc(M, eta, H, (const ConstCMD<Type>&)*this);  Type MxH[3]; vprod(M, H, MxH);  			
			// WMSG(coeff.mu_p, coeff.nu_p, coeff.THETA[0],coeff.THETA[1],coeff.THETA[2], coeff.Mabs, coeff.Psi, coeff.U, coeff.Q, coeff.MnK);
			for(int i=0; i<3; i++) dM[i] = -MxH[i] + 2*K*(coeff.PHI[i] + alpha*coeff.THETA[i]) + alpha*(coeff.XiH[i] + 2*(n_b*J*coeff.U - T)*M[i]);
			deta = 4*alpha*(sprod(M, H)*coeff.U + K*coeff.Psi + J*coeff.Q - T*eta);   // проправка с dMdt в MU?
			if(coeff.zeta>coeff.zeta_c){
				Type C = CdMdt*(coeff.zeta-coeff.zeta_c);
				deta += 4*alpha*C*sprod(dM, H);
				for(int i=0; i<3; i++) dM[i] += 2*alpha*n_b*J*C*dM[i];
			}
		}
		//----------------------------------------------------------------------		
	};
	//--------------------------------------------------------------------------
	//  shot noise
	//--------------------------------------------------------------------------
	template <typename T, typename R> T rand_M(T M[3], T sigma, R& N01gen){  // sigma = sqrt(2*alpha*gamma*T*dt/N);
		T M2 = M[0]*M[0]; for(int i=1; i<3; i++) M2 += M[i]*M[i];
		T rnd[3];  for(int i=0; i<3; i++) rnd[i] = N01gen();
		T Mabs = sqrtf(M2), Mabs2;  
		if(Mabs<sigma){
			const T sqrt23 = sqrtf(2.f/3); M2 = 0;
			for(int i=0; i<3; i++){ M[i] += rnd[i]*sigma*sqrt23; M2 += M[i]*M[i]; }
			Mabs2 = sqrtf(M2);
		} else {
			T scale = Mabs; 
			// T scale = Mabs + sigma*sqrtf(.666666f*(1-M2)*(1+.488f*M2))*sprod(M, rnd); 
			T Mxrnd[3];  vprod(M, rnd, Mxrnd);	M2 = 0;
			for(int i=0; i<3; i++){ M[i] += Mxrnd[i]*sigma;  M2 += M[i]*M[i]; }
			Mabs2 = sqrtf(M2);  scale /= Mabs2; if(fabsf(scale*Mabs)>1){ scale = 1/Mabs2; Mabs2 = 1; }
			for(int i=0; i<3; i++) M[i] *= scale; 
		}
		return Mabs2-Mabs;
	}
	//--------------------------------------------------------------------------
	template <typename T, typename R> void rand_M_eta(T M[3], T &eta, T sigma, R& N01gen){
		rand_M(M, sigma, N01gen); 
		// eta += rand_M(M, sigma, N01gen) + sigma*(1-eta)*N01gen(); 
		// eta += dM + sigma*(eta<1? sqrtf(1-eta): 0)*N01gen();
		// eta += -dM + sigma*(1-eta)*N01gen();
		// T s = 0.97f-1.77*eta+5.39f*eta*eta-4.35*eta*eta*eta;
		// eta +=  sigma*s*N01gen();
		// T eta2 = eta*eta;
		// eta += sqrtf(sigma*.25f*M_PI*(1+(-.444f+(.253f-.612f*eta2)*eta2)*eta2))*N01gen();
		// if(eta<0) eta = 0;
		// else if(eta>1) eta = 1;  // ??? 1-eps ???
	}
	//--------------------------------------------------------------------------
}

#endif // CMDLIB_CMD_COEFF_HPP
