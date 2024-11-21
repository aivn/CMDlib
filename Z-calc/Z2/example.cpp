#include <fstream>
#include <iostream>

#include "Z2.hpp"

int main() {
    double min  = 0.01;
    double max  = 0.99;
    double N    = 100;
    double step = (max - min) / N;

    Z2 m_Z2 = Z2();

    std::ofstream out;
    out.open("Z2_coeffs.dat");
    out << "#:m eta h l steps" << '\n';

    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            double m = min + j * step, eta = min + i * step;
            
            m_Z2.calc(m, eta);
            if (m*m < eta)
                out << m << ' ' << eta << ' ' << m_Z2.h << ' ' << m_Z2.lambda << ' ' << m_Z2.steps << '\n';
        }
        out << '\n';
    }
    out.close();
}
