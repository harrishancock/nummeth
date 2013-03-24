#include "polynomial.hpp"

#include <cstdio>

int main () {
    polynomial p { { 2.5, 0, 12 } };
    printf("p(1) == %f\n", p(1));
    printf("p(3) == %f\n", p(3));
    printf("p(8) == %f\n", p(8));
    polynomial dp { p.derivative() };
    printf("dp(1) == %f\n", dp(1));

    polynomial p0 { { 1, -12, 0, -42 } };
    polynomial p1 { { 1, -3 } };
    polynomial p2 { { 1, 1, -3 } };

    polynomial p3 { { 6, 5, 0, -7 } };
    polynomial p4 { { 3, -2, -1 } };

    {
        auto qr = p0 / p1;
        const polynomial& q = qr.first;
        const polynomial& r = qr.second;
        std::cout << p0 << " / " << p1 << " == (" << q << ", " << r << ")\n";
    }

    {
        auto qr = p0 / p2;
        const polynomial& q = qr.first;
        const polynomial& r = qr.second;
        std::cout << p0 << " / " << p2 << " == (" << q << ", " << r << ")\n";
    }

    {
        auto qr = p3 / p4;
        const polynomial& q = qr.first;
        const polynomial& r = qr.second;
        std::cout << p3 << " / " << p4 << " == (" << q << ", " << r << ")\n";
    }

    return 0;
}
