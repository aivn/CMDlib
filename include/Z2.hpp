/**
 * @file Z2.hpp
 * @author Andrei V. Lukianov <andrey.luk10@gmail.com>
 * @brief Расчет двухчастичной статистической суммы.
 * 
 * Симметричный случай: случай, при котором hi = hj = h.
 * 'norm' в симметричном случае означает умножение на exp(-2h-l).
 * 
 * Несимметричный случай: случай, при котором hi != hj, но hi и hj коллинеарны.
 * 'norm' в несимметричном случае означает умножение на exp(- hi-hj-l).
 * 
 * @copyright Copyright (c) 2025
 * Licensed under the Apache License, Version 2.0
 */
#ifndef CMDLIB_Z2_HPP
#define CMDLIB_Z2_HPP

// =======================================================================================

inline double L(double p) { return fabs(p) > 1e-1 ? 1 / tanh(p) - 1 / p : p / 3; }
inline double dLdp(double p) {
    double shp = sinh(p);
    return fabs(p) > 1e-6 ? 1 / (p * p) - 1 / (shp * shp) : 1 / 3. - p * p / 15;
}
inline double invL(double M) {
    if (fabsf(M) < 1e-6) return M * 3;
    double p = 1.f;
    for (int i = 0; i < 10; i++) {
        double Lp = L(p);
        if (p > 400 || fabs(Lp - M) < 1e-6) break;
        p = p - (Lp - M) / dLdp(p);
        if (fabs(p) < 1e-6) return p;
    }
    return p;
}

// =======================================================================================

double cF2_norm_symmetrical(double h, double l);
double cF2_symmetrical(double h, double l);

double E_norm_symmetrical(double h, double l, bool arg_s1, bool arg_s2);
double E_symmetrical(double h, double l, bool arg_s1, bool arg_s2);

/**
 * @brief Структура для решение прямой и обратной задач для 
 * двухчастичной статистической суммы в симметричном случае.
 * Ограничения: eta >= m*m, l >= 0.
 */
struct Z2_symmetrical {
private:
    double _h, _l, sqrt_l, _sqrt_l;
    double _exp_2h, _exp_2l;
    double cF2_norm, E12_norm, E1_norm, E_norm;
    // double m_llb, m_h_llb, eta_prmg;
    double m_h;

    void calc_base();

    void calc_m();
    void calc_eta();
    void calc_upsilon();
    void calc_eta2();
    void calc_mh2();
    void calc_psi0();
    void calc_Z2_norm();

public:
    double h, l;
    double m, eta;
    double m_eta, upsilon, eta2;
    double mh2, psi0;
    double Z2_norm;

    /**
     * @brief Рассчитывает прямую задачу: h, l -> moments.
     * Значение h, l задаются через поля.
     * Ограничения: l >= 0.
     * 
     */
    void calc_from_coeffs();

    /**
     * @brief Решает обратную задачу: m, eta -> h, l -> moments.
     * Значение m, eta задаются через поля.
     * Ограничения: eta >= m*m, eta > 0. 
     * Значение h, l хранятся как поля.
     * После нахождения h, l рассчитывает все моменты (m, eta пересчитываются).
     * 
     * @param h0 Начальное значение h.
     * @param l0 Начальное значение l.
     * @param eps Невязка.
     * @param max_steps Максимально допустимое количество шагов.
     * @return int Количество шагов метода Ньютона.
     */
    int calc_from_moments(double h0, double l0, double eps = 1e-8, int max_steps = 64);
};


double cF2_norm_not_symmetrical(double hi, double hj, double l);
double cF2_not_symmetrical(double hi, double hj, double l);

double E_norm_not_symmetrical(double hi, double hj, double l, bool arg_s1, bool arg_s2);
double E_not_symmetrical(double hi, double hj, double l, bool arg_s1, bool arg_s2);

#endif  // CMDLIB_Z2_HPP
