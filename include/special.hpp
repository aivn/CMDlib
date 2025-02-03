/**
 * @file special.hpp
 * @author Andrei V. Lukianov <andrey.luk10@gmail.com>
 * @brief Файл с реализациями спецфункций.
 * 
 * @copyright Copyright (c) 2025
 * Licensed under the Apache License, Version 2.0
 */
#ifndef CMDLIB_SPECIAL_HPP
#define CMDLIB_SPECIAL_HPP

inline int sign(float val) {
    return (0 < val) - (val < 0);
}

/**
 * @brief Реализация функции Доусона. Абсолютная точность около 2e-7.
 */
float dawson(float x);

/**
 * @brief Реализация scaled complementary error function для x >= 0 с относительной ошибкой < 1.2e-7.
 */
float erfcx(float x);

/**
 * @brief Реализация функции F3. Абсолютная точность около 1e-8.
 */
float F3(float x, float a, float b);

#endif  // CMDLIB_SPECIAL_HPP

