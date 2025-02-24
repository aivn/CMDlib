/**
 * @file special.hpp
 * @author Andrei V. Lukianov <andrey.luk10@gmail.com>
 * @brief Реализации необходимых спецфункций.
 * 
 * @copyright Copyright (c) 2025
 * Licensed under the Apache License, Version 2.0
 */
#ifndef CMDLIB_SPECIAL_HPP
#define CMDLIB_SPECIAL_HPP

/**
 * @brief Реализация функции Доусона. Абсолютная точность около 2e-7.
 */
double dawson(double x);

/**
 * @brief Реализация scaled complementary error function для x >= 0.
 */
double erfcx(double x);

/**
 * @brief Реализация функции F3. Абсолютная точность около 1e-7.
 */
double F3(double x, double a, double b);

#endif  // CMDLIB_SPECIAL_HPP

