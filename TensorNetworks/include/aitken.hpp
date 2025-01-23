#ifndef AITKEN_HPP
#define AITKEN_HPP

#include "math.h"
#include <vector>

/**
 * @brief Интерполяцию в точку 0.
 */
template <typename TF, typename TA>
TF interpolation(const TF* G, const TA* X, int N, TA x){
    std::vector<TF> GG(N);
    for(int i = 0; i < N-1;i++){
        TA x1 = X[i]; TA x2 = X[i+1];
        GG[i] = (G[i+1]*(x-x1)+G[i]*(x2-x))/(x2-x1);
	}
	for(int step=2; step < N; step++){
		for(int i = 0; i < N-step;i++){
			TA x1 = X[i];
			TA x2 = X[i+step];
			GG[i] = (GG[i+1]*(x-x1)+GG[i]*(x2-x))/(x2-x1);
		}
	}
	return GG[0];
}

/**
 * @brief Функция, реализующая итерационное вычисление интеграла по Эйткену.
 */
template<typename TF, typename TA, typename TR>
double aitken(TF func, TA x_0, TA x, TR eps, int N = 10){
    TA a = x_0;
	TA b = x;
	TA h=(b-a)/N;
	std::vector<TA> H;
	std::vector<TR> I_;
	H.push_back(h*h);
	TR I_last=0.5*(func(a) + func(b));
	for(int i=1;i<N;i++){
		I_last += func(a+i*h);
	}
	I_last *= h;
	I_.push_back(I_last);
	
	TR dif;
	int count = 0;
	do{
		h = 0.5*h;
		N = N*2;
		H.push_back(h*h);
		TR I_cur = 0.;
		for(int i=0;i<N;i+=2){
			I_cur += func(a+h+i*h);
		}
		I_cur *= h;
		I_cur += I_[I_.size()-1]*0.5;
		I_.push_back(I_cur);
		
		I_cur = interpolation(&I_[0], &H[0], int(I_.size()), (TA)0);
		dif = std::abs(I_last-I_cur);
		I_last = I_cur;
		count++;
	} while(dif > eps);

	return I_last;
}

#endif //AITKEN_HPP