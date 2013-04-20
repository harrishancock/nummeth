#include "polynomial.hpp"

#include <cstdio>

int main () {
    polynomial p { { 12, 0, 2.5 } };
    printf("p(1) == %f\n", p(1));
    printf("p(1.1) == %f\n", p(1.1));
    printf("p(1.2) == %f\n", p(1.2));
    printf("p(1.3) == %f\n", p(1.3));
    printf("p(3) == %f\n", p(3));
    printf("p(8) == %f\n", p(8));
    polynomial dp { p.derivative() };
    printf("dp(1) == %f\n", dp(1));

    polynomial p0 { { -42, 0, -12, 1 } };


    polynomial p1 { { -3, 1 } };
    polynomial p2 { { -3, 1, 1 } };

    polynomial p3 { { -7, 0, 5, 6 } };
    polynomial p4 { { -1, -2, 3 } };

    std::cout << "p0 == " << p0 << '\n';
    std::cout << "p1 == " << p1 << '\n';
    std::cout << "p2 == " << p2 << '\n';
    std::cout << "p3 == " << p3 << '\n';
    std::cout << "3*p0 == " << 3*p0 << '\n';
    std::cout << "p0 + p1 == " << p0 + p1 << '\n';
    std::cout << "p0 + p2 == " << p0 + p2 << '\n';
    std::cout << "p0 + p3 == " << p0 + p3 << '\n';
    std::cout << "p0 - p1 == " << p0 - p1 << '\n';
    std::cout << "p0 - p2 == " << p0 - p2 << '\n';
    std::cout << "p0 - p3 == " << p0 - p3 << '\n';
    std::cout << "p0 * p1 == " << p0 << " * " << p1 << " == " << p0 * p1 << '\n';
    std::cout << "p0 * p2 == " << p0 << " * " << p2 << " == " << p0 * p2 << '\n';
    std::cout << "p0 * p3 == " << p0 << " * " << p3 << " == " << p0 * p3 << '\n';

    {
        auto qr = p0.quo_rem(p1);
        const polynomial& q = qr.first;
        const polynomial& r = qr.second;
        std::cout << p0 << " / " << p1 << " == (" << q << ", " << r << ")\n";
    }

    {
        auto qr = p0.quo_rem(p2);
        const polynomial& q = qr.first;
        const polynomial& r = qr.second;
        std::cout << p0 << " / " << p2 << " == (" << q << ", " << r << ")\n";
    }

    {
        auto qr = p3.quo_rem(p4);
        const polynomial& q = qr.first;
        const polynomial& r = qr.second;
        std::cout << p3 << " / " << p4 << " == (" << q << ", " << r << ")\n";
    }

    polynomial sturm_p { { -1, -1, 0, 1, 1 } };
    std::cout << "Sturm Chain of (" << sturm_p << ")\n";

    auto sturm_c = sturm_chain(sturm_p);
    for (auto p : sturm_c) {
        std::cout << p << '\n';
    }

    {
        printf("p[i] at -inf: ");
        int sigma = 0;
        bool last_sign = sturm_c[0][sturm_c[0].degree()] < 0;
        last_sign = sturm_c[0].degree() & 1 ? !last_sign : last_sign;
        for (auto p : sturm_c) {
            bool sign = p[p.degree()] < 0;
            sign = p.degree() & 1 ? !sign : sign;
            sigma += last_sign != sign;
            last_sign = sign;
            printf("%c ", sign ? '-' : '+');
        }
        printf("(%d)\n", sigma);
    }

    {
        printf("p[i] at inf: ");
        int sigma = 0;
        bool last_sign = sturm_c[0][sturm_c[0].degree()] < 0;
        for (auto p : sturm_c) {
            bool sign = p[p.degree()] < 0;
            sigma += last_sign != sign;
            last_sign = sign;
            printf("%c ", sign ? '-' : '+');
        }
        printf("(%d)\n", sigma);
    }

    {
        printf("Legendre polynomials:\n");
        polynomial x { { 0, 1 } };
        std::vector<polynomial> legendre (10, polynomial { { 0.0 } });
        legendre[0] = polynomial { { 1 } };
        legendre[1] = x;
        for (int i = 2; i < legendre.size(); i++) {
            legendre[i] = ((2*i - 1)*x*legendre[i-1] - (i - 1)*legendre[i-2]) * (1.0/i);
        }
        for (auto p : legendre) {
            std::cout << p << '\n';
        }
    }

    return 0;
}
