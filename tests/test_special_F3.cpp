#include <cassert>
#include <iostream>

#include <math.h>

#include "../include/special.hpp"

int main() {
    assert(F3(0, 0, 0) == 0.);
    assert(fabs(F3(1, .5, 1) - 0.25897852311227765) < 1e-7);
    assert(fabs(F3(1, 1.5, 1) - 0.18357150587846913) < 1e-7);
    assert(fabs(F3(1, -9.90399 , 1) + 0.06014129314320898) < 1e-7);
    assert(fabs(F3(1, -10 + 56 * 20. / 99., 2) - 0.10927116426850177) < 2e-7);

    int N       = 1000;
    double left = -10, right = 10;
    double h = (right - left) / (N - 1);
    for (int i = 0; i < N; i++){
        double x = 1;
        double a = .5;
        double b = left + i*h;
        std::cout << b << ' ' << F3(x, a, b) << '\n';
    }

    return 0;
}
