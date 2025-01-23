#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <vector>
#include <ostream>
#include <math.h>

template<typename T>
class Polynomial{
    int n;
    std::vector<T> coeffs;
public:
    int get_n() const { return n; }
    Polynomial(int n_ = 0){
        n = n_;
        coeffs.resize(n+1);
    }
    template <typename other_Type> Polynomial<T>& operator = (const Polynomial<other_Type> &p);
    T calc(T x) const{
        T res = 0;
        for(int i = 0; i < n+1; i++){ res += coeffs[i]*pow(x, i); }
        return res;
    }

    T & operator[](int i){ return coeffs[i]; }
    T operator[] (int i) const{ return coeffs[i]; }

    Polynomial<T> diff(int k = 1) const{
        if(k == 0){
            return *this;
        }else{
            Polynomial<T> res(n-k);
            for(int i = 0; i < n-k+1; i++){
                T & v = res[i];
                v = coeffs[i+k];
                for(int j = 0; j < k; j++){
                    v *= i+k-j;
                }
            }
            return res;
        }
    }
};

template <typename T1, typename T2> Polynomial<decltype(T1()+T2())> operator + (const Polynomial<T1> & p1, const Polynomial<T2> & p2){
    int n1 = p1.get_n(); int n2 = p2.get_n();
    int max_n = std::max(n1, n2);
    int min_n = std::min(n1, n2);
    int max_pol = max_n==n1 ? 0:1;
    Polynomial<decltype(T1()+T2())> res(max_n);
    for(int i = 0; i < min_n+1; i++){
        res[i] = p1[i]+p2[i];
    }
    
    for(int i = min_n+1; i < max_n+1; i++){
        res[i] = max_pol==0? p1[i]:p2[i];
    }
    return res;
}

template <typename T1, typename T2> Polynomial<decltype(T1()-T2())> operator - (const Polynomial<T1> & p1, const Polynomial<T2> & p2){
    return p1+(-1)*p2;
}

template <typename T1, typename T2> Polynomial<decltype(T1()*T2())> operator * (const Polynomial<T1> & p1, const Polynomial<T2> & p2){
    int n1 = p1.get_n(); int n2 = p2.get_n(); int n = n1*n2;
    Polynomial<decltype(T1()*T2())> res(n+1);
    
    for(int i = 0; i < n1+1; i++){
        for(int j = 0; j < n2+1; j++){
            res[i+j] += p1[i]*p2[j];
        }
    }
    return res;
}

template <typename Type>
template <typename other_Type> Polynomial<Type>& Polynomial<Type>::operator = (const Polynomial<other_Type> &p){
    n = p.get_n();
    coeffs.resize(n+1);
    for(int i = 0; i < n+1; i++){
        coeffs[i] = p[i];
    }
    return *this;
}

template <typename T>
std::ostream & operator << (std::ostream & S, const Polynomial<T> & p){
    S << "[";
    int n = p.get_n(); 
    for(int i = 0; i < n; i++){
        S << i << ":" << p[i] << ", ";
    }
    S << n << ":" << p[n] << "]";
    return S;
}

template <typename T1, typename T2>
Polynomial<decltype(T1()*T2())>  operator * (const Polynomial<T1> & p, T2 a){
    Polynomial<decltype(T1()*T2())> res = p;
    int n = p.get_n();
    for(int i = 0; i < n+1; i++){
        res[i] *= a;
    } 
    return res;
}

template <typename T1, typename T2>
Polynomial<decltype(T1()*T2())>  operator * (T2 a, const Polynomial<T1> & p){
    Polynomial<decltype(T1()*T2())> res = p;
    int n = p.get_n();
    for(int i = 0; i < n+1; i++){
        res[i] *= a;
    } 
    return res;
}

#endif //POLYNOMIAL_HPP