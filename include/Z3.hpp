/**
 * @file Z2.hpp
 * @author Andrei V. Lukianov <andrey.luk10@gmail.com>
 * @brief Расчет трехчастичной статистической суммы.
 * 
 * Симметричный случай: случай, при котором hi = hk и l_ij = l_jk.
 * 'norm' в симметричном случае означает умножение на exp(-h_j-l_ik-2*(h_i+l_ij)).
 * 
 * @copyright Copyright (c) 2025
 * Licensed under the Apache License, Version 2.0
 */
#ifndef CMDLIB_Z3_HPP
#define CMDLIB_Z3_HPP

double cF3_norm_symmetrical(double h_i, double h_j, double l_ij, double l_ik);
double cF3_symmetrical(double h_i, double h_j, double l_ij, double l_ik);

double E_norm_symmetrical(double h_i, double h_j, double l_ij, double l_ik, bool arg_s1, bool arg_s2, bool arg_s3);
double E_symmetrical(double h_i, double h_j, double l_ij, double l_ik, bool arg_s1, bool arg_s2, bool arg_s3);

double FTilde_norm_symmetrical(double h_i, double h_j, double l_ij, double l_ik, bool arg_s1, bool arg_s2, bool arg_s3);
double FTilde_symmetrical(double h_i, double h_j, double l_ij, double l_ik, bool arg_s1, bool arg_s2, bool arg_s3);

double F_norm_symmetrical(double h_i, double h_j, double l_ij, double l_ik, bool arg_s1, bool arg_s2, bool arg_s3);
double F_symmetrical(double h_i, double h_j, double l_ij, double l_ik, bool arg_s1, bool arg_s2, bool arg_s3);

/**
 * @brief Структура для решение прямой и обратной задач для 
 * трехчастичной статистической суммы в симметричном случае.
 * Ограничения для прямой задачи: l_ij > 0, l_ik > 0.
 * Ограничения для обратной задачи: mi = mj = m, eta_ij > m*m, eta_ik > m*m, eta_ik > eta_ij*eta_ij
 */
// struct Z3_symmetrical {
// private:
//     double _h, _l, sqrt_l, _sqrt_l;
//     double _exp_2h, _exp_2l;
//     double cF2_norm, E12_norm, E1_norm, E_norm;
//     double m_llb, m_h_llb, eta_prmg;
//     double m_h;

//     void calc_base();

//     void calc_mi();
//     void calc_mj();
//     void calc_eta_ij();
//     void calc_eta_ik();
//     void calc_mj_eta_ij();
//     void calc_mjh2();
//     void calc_Qstar();

// public:
//     double h_i, h_j, l_ij, l_ik;
//     double mi, mj, m, eta_ij, eta_ik;
//     double mj_eta_ij, mjh2, Qstar;
//     double Z2_norm;

//     /**
//      * @brief Рассчитывает прямую задачу: coeffs -> moments.
//      * Значение коэффициентов задаются через поля.
//      * Ограничения: l_ij > 0, l_ik > 0.
//      * 
//      */
//     void calc_from_coeffs();

//     /**
//      * @brief Решает обратную задачу: m, eta_ij, eta_ij -> h, l -> moments.
//      * Значение m (не mi и mj (!)), eta_ij, eta_ij задаются через поля.
//      * Ограничения: mi = mj = m, eta_ij > m*m, eta_ik > m*m, eta_ik > eta_ij*eta_ij.
//      * Значение коэффициентов хранятся как поля.
//      * После нахождения коэффициентов рассчитывает все моменты (mi, mj, eta_ij, eta_ij пересчитываются).
//      * 
//      * @param hi Начальное значение hi.
//      * @param hj Начальное значение hj.
//      * @param lij Начальное значение lij.
//      * @param lik Начальное значение lik.
//      * @param eps Невязка.
//      * @param max_steps Максимально допустимое количество шагов.
//      * @return int Количество шагов метода Ньютона.
//      */
//     int calc_from_moments(double hi0, double hj0, double lij0, double lik0, double eps = 1e-8, int max_steps = 64);
// };

#endif  // CMDLIB_Z3_HPP
