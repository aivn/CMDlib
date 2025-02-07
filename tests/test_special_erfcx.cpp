#include <cassert>
#include <iostream>

#include <math.h>

#include "../include/special.hpp"

int main() {
    assert(erfcx(0) == 1.);
    assert(fabs(erfcx(1) - 0.427583576155807) < 1.2e-7);

    int N = 100;
    double left = 0, right = 10;
    double h = (right-left)/(N-1);
    for (int i = 0; i < N; i++){
        double x = left + i*h;
        std::cout << x << ' ' << erfcx(x) << '\n';
    }

    return 0;
}
