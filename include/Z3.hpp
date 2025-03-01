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
 * Ограничения: l_ij > 0, l_ik > 0.
 */
struct Z3_symmetrical {
private:
    double kappa, _kappa;
    double sqrt_h_i, sqrt_h_j, sqrt_l_ij, sqrt_l_ik;
    double _exp_2h_i, _exp_2h_j, _exp_2l_ij, _exp_2l_ik;
    double cF3, _cF3;
    double E123, E12, E1, E13, E3, E;
    double FTilde3, FTilde123, FTilde12;
    double F13, F, F1;

    void calc_base();

    void calc_mi();
    void calc_mj();
    void calc_eta_ij();
    void calc_eta_ik();

    void calc_mj_eta_ij();
    void calc_mi_eta_ik();
    void calc_mj_eta_ik();
    void calc_eta_ik2();

    void calc_mjh2();
    void calc_Q_star();

    void calc_Z3_norm();

public:
    double h_i, h_j, l_ij, l_ik;
    double mi, mj, eta_ij, eta_ik;
    double mj_eta_ij, mi_eta_ik, mj_eta_ik, eta_ik2;
    double mjh2, Q_star;
    double Z3_norm;

    /**
     * @brief Рассчитывает прямую задачу: coeffs -> moments.
     * Значение коэффициентов задаются через поля.
     * Ограничения: l_ij > 0, l_ik > 0.
     */
    void calc_from_coeffs();

    /**
     * @brief Решает обратную задачу: mi, mj, eta_ij, eta_ij -> h, l -> moments.
     * Значение mi, mj, eta_ij, eta_ij задаются через поля.
     * Ограничения:l_ij > 0, l_ik > 0.
     * Значение коэффициентов хранятся как поля.
     * После нахождения коэффициентов рассчитывает все моменты (mi, mj, eta_ij, eta_ij пересчитываются).
     * 
     * @param hi Начальное значение hi.
     * @param hj Начальное значение hj.
     * @param lij Начальное значение lij.
     * @param lik Начальное значение lik.
     * @param eps Невязка.
     * @param max_steps Максимально допустимое количество шагов.
     * @return int Количество шагов метода Ньютона.
     */
    int calc_from_moments(double hi0, double hj0, double lij0, double lik0, double eps = 1e-8, int max_steps = 64);
};

#endif  // CMDLIB_Z3_HPP
