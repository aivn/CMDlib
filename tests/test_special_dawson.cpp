#include <cassert>
#include <iostream>

#include <math.h>

#include "../include/special.hpp"

int main() {
    assert(dawson(0) == 0.);
    assert(fabs(dawson(1) - 0.5380795069127684) < 2e-7);
    assert(dawson(1) == -dawson(-1));

    int N = 100;
    double left = -10, right = 10;
    double h = (right-left)/(N-1);
    for (int i = 0; i < N; i++){
        double x = left + i*h;
        std::cout << x << ' ' << dawson(x) << '\n';
    }

    return 0;
}
