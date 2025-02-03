#include "../include/special.hpp"

#include <algorithm>
#include <math.h>

float dawson(float x) {
    static const int NMAX = 6;
    static const float h = .4, A1 = 2.0 / 3.0, A2 = 0.4, A3 = 2.0 / 7.0;

    // exp(-((2*i+1)*h)**2), 0 <= i < NMAX
    static const float coeffs[6] = {
        0.8521437889662113, 0.23692775868212165, 0.01831563888873418,
        0.0003936690406550776, 2.352575200009771e-06, 3.90893843426485e-09
    };

    if (abs(x) < 0.2) {  // Разложение в ряд.
        float x2 = x * x;
        return x * (1.0 - A1 * x2 * (1.0 - A2 * x2 * (1.0 - A3 * x2)));
    } else {  // Sampling theorem.
        float xx = fabs(x);
        int n0   = 2 * int(0.5 * xx / h + 0.5);
        float xp = xx - n0 * h;

        float e1 = exp(2.0 * xp * h);
        float e2 = e1 * e1;
        float d1 = n0 + 1;
        float d2 = d1 - 2.0;

        float result = 0.0;
        for (int i = 0; i < NMAX; i++) {
            result += coeffs[i] * (e1 / d1 + 1.0 / (d2 * e1));
            d1 = d1 + 2.0;
            d2 = d2 - 2.0;
            e1 = e1 * e2;
        }

        return 0.5641895835 * exp(-xp * xp) * sign(x) * result;
    }
}

float erfcx(float x) {
    if (x < 0.) throw("erfcx requires nonnegative argument");

    double t = 2. / (2. + x);

    double result = t * exp(-1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return result;
}

/**
 * @brief Антисимметричный вариант функции erfcx. Имеет разрыв в нуле.
 */
float antisymmetric_erfcx(float x) {
    if (x >= 0) return erfcx(x);
    else return -erfcx(-x);
}

const float __F3_h = .4;

/**
 * @brief Функция, рассчитывающая значения сумм S_1 и S_2 (см. текст).
 * 
 * @param func Для c > 0 -- это dawson, для c < 0 -- это antisymmetric_erfcx. 
 * Принимает на вход значение аргумента n.
 * 
 * @param exp_arg0 Часть аргумента экспоненты, не зависящий от n.
 * @param exp_arg_delta Часть аргумента экспоненты, линейно зависящая от n.
 */
template <typename TFunc>
float __F3_calc_sum(int n0, TFunc func, float exp_arg0, float exp_arg_delta) {
    static const int NMAX = 6;

    // exp(-((2*i+1)*h)**2), 0 <= i < NMAX
    static const float coeffs[6] = {
        0.8521437889662113, 0.23692775868212165, 0.01831563888873418,
        0.0003936690406550776, 2.352575200009771e-06, 3.90893843426485e-09
    };

    int n1 = n0 + 1;
    int n2 = n0 - 1;

    float e1        = exp(exp_arg0 + .5 * exp_arg_delta);
    float e2        = exp(exp_arg0 - .5 * exp_arg_delta);
    float delta_exp = exp(exp_arg_delta);

    float result = 0;
    for (int i = 0; i < NMAX; i++) {
        result += func(n1) * e1 * coeffs[i] / n1;
        result += func(n2) * e2 * coeffs[i] / n2;

        n1 += 2;
        n2 -= 2;
        e1 *= delta_exp;
        e2 /= delta_exp;
    }

    return result;
}

/**
 * @brief Функция, рассчитывающая значения суммы S_1 (см. текст).
 * @param func см __F3_calc_sum.
 */
template <typename TFunc>
float __F3_calc_S1(TFunc func, float x, float a, float b, float c) {
    float _sqrt_c = 1. / sqrt(fabs(c));

    float h  = __F3_h;
    float xx = (a * x + b) / h;
    int n0   = 2 * int(0.5 * xx + sign(xx) * 0.5);

    float arg0_exp      = 2 * b * h * n0 - h * h * n0 * n0 + 2 * a * h * n0 * x - 2 * a * b * x + (c - 1) * x * x - b * b;
    float delta_arg_exp = 4 * a * h * x + 4 * b * h - 4 * h * h * n0;

    float func_arg_linear = a * h * _sqrt_c;
    float func_arg_const  = (-a * b + c * x) * _sqrt_c;

    float zero_n      = (c * x + a * b) / (a * h);
    auto applied_func = [zero_n, func_arg_linear, func_arg_const, func](int n) {
        if (n == zero_n) return func(0.);
        return func(func_arg_linear * n + func_arg_const);
    };

    return __F3_calc_sum(n0, applied_func, arg0_exp, delta_arg_exp);
}

/**
 * @brief Функция, рассчитывающая значения суммы S_2 (см. текст).
 * @param func см __F3_calc_sum.
 */
template <typename TFunc>
float __F3_calc_S2(TFunc func, float x, float a, float b, float c) {
    float _sqrt_c = 1. / sqrt(fabs(c));

    float h  = __F3_h;
    float xx = b / h;
    int n0   = 2 * int(0.5 * xx + sign(xx) * 0.5);

    float arg0_exp      = 2 * b * h * n0 - h * h * n0 * n0 - b * b - x * x;
    float delta_arg_exp = 4 * b * h - 4 * h * h * n0;

    float func_arg_linear = a * h * _sqrt_c;
    float func_arg_const  = -a * b * _sqrt_c;

    float zero_n      = b / h;
    auto applied_func = [zero_n, func_arg_linear, func_arg_const, func](int n) {
        if (n == zero_n) return func(0.);
        return func(func_arg_linear * n + func_arg_const);
    };

    return __F3_calc_sum(n0, applied_func, arg0_exp, delta_arg_exp);
}

float F3(float x, float a, float b) {
    float h = __F3_h;
    float c = 1 - a * a;

    if (c == 1)
        return dawson(x) * dawson(b);

    if (c > 0) {
        float _sqrt_c = 1 / sqrt(c);

        float result = __F3_calc_S1([](float x) { return dawson(x); }, x, a, b, c);
        result -= __F3_calc_S2([](float x) { return dawson(x); }, x, a, b, c);

        return result / sqrt(M_PI) * _sqrt_c;
    }

    if (c < 0) {
        float _sqrt_c = 1 / sqrt(-c);

        float result = __F3_calc_S1([](float x) { return antisymmetric_erfcx(x); }, x, a, b, c);
        result -= __F3_calc_S2([](float x) { return antisymmetric_erfcx(x); }, x, a, b, c);

        float factor = -2 * sign(a);
        float left = (a * b - c * x) / (a * h), right = b / h;

        if (left != right) {
            if (right < left) {
                factor *= -1;
                std::swap(left, right);
            }

            int n_prob  = (left + right) / 2;
            float start = left, end = right;
            if (a * a * (b / c - h * n_prob / c) + b - h * n_prob > 0)
                std::swap(start, end);

            int step = 2 * sign(end - start);
            int n    = floor(start);
            n        = n - 1 + abs(n) % 2;

            float distance = fabs(end - start);
            if (fabs(n - end) > distance)
                n += step;

            while (true) {
                float arg_exp = (-a * a * b * b + 2 * a * a * b * h * n - a * a * h * h * n * n) / c - b * b + 2 * b * h * n - h * h * n * n - x * x;
                if (fabs(n - start) > distance || fabs(n - end) > distance) break;
                if (n == left)
                    if (a < 0) result -= factor * exp(arg_exp) / n;
                if (n == right) {
                    if (a < 0) result += factor * exp(arg_exp) / n;
                    n += step;
                    continue;
                }

                float val_exp = exp(arg_exp);
                result += factor * val_exp / n;
                n += step;

                if (fabs(factor * val_exp / n) - 1e-7 * fabs(result) < 0) break;
            }
        }

        return result * .5 * _sqrt_c;
    }
    // c == 0
    if (b == 0) return x * dawson(x) - .5 * (1 - exp(-x * x));

    // b != 0
    float result = dawson(b) * exp(-x * x);
    result += dawson(x);
    result -= dawson(b + x);
    return result * .5 / b;
}
