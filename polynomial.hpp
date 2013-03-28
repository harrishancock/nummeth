#ifndef NUMMETH_POLYNOMIAL_HPP
#define NUMMETH_POLYNOMIAL_HPP

#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/combine.hpp>

#include <vector>
#include <iostream>

class polynomial {
public:
    polynomial () : coefficients{0} { }

    explicit polynomial (std::vector<double> cs) : coefficients(cs) { }

    double operator() (const double x) const {
        double p = 0;
        for (auto c : coefficients | boost::adaptors::reversed) {
            p = p * x + c;
        }
        return p;
    }

    double operator[] (const size_t i) const {
        return i > degree() ? 0 : coefficients[i];
    }

    size_t degree () const {
        return coefficients.size() - 1;
    }

    polynomial derivative () const {
        if (!degree()) {
            return polynomial { { 0 } };
        }

        std::vector<double> dpdx { coefficients.cbegin() + 1, coefficients.cend() };

        double power = 1.0;
        for (auto& c : dpdx) {
            c *= power++;
        }
        
        return polynomial { dpdx };
    }

    template <typename DIter, typename QIter>
    static auto ip (DIter dbegin, DIter dend, QIter qbegin, QIter qend) -> decltype(*dbegin * *qbegin) {
        decltype(*dbegin * *qbegin) value = 0;
        while (dbegin != dend && qbegin != qend) {
            value += *dbegin * *qbegin;
            ++dbegin;
            ++qbegin;
        }
        return value;
    }

    std::pair<polynomial, polynomial> quo_rem (polynomial divisor) const {
        /* Degrees of the remainder and quotient. */
        int kd = divisor.degree();
        int kr = kd - 1;
        int kq = degree() - kd;

        std::vector<double> remainder (kr + 1);
        std::vector<double> quotient (kq + 1);

        auto dbegin = divisor.coefficients.crbegin() + 1;
        auto dend = divisor.coefficients.crend();
        auto qbegin = quotient.end();
        auto qend = quotient.end();

        for (int i = kq; 0 <= i; i--) {
            quotient[i] = (coefficients[kd+i] - ip(dbegin, dend, qbegin, qend)) / divisor[kd];
            --qbegin;
        }

        for (int i = kr; 0 <= i; i--) {
            remainder[i] = coefficients[i] - ip(dbegin, dend, qbegin, qend);
            ++dbegin;
        }

        return { polynomial { quotient }, polynomial { remainder } };
    }

    const polynomial operator/ (polynomial divisor) const {
        return quo_rem(divisor).first;
    }

    const polynomial operator% (polynomial divisor) const {
        return quo_rem(divisor).second;
    }

    const polynomial operator+ (polynomial operand) const {
        std::vector<double> result (std::max(degree(), operand.degree()) + 1);

        for (int i = 0; i < result.size(); i++) {
            double x = i <= degree() ? coefficients[i] : 0.0;
            double y = i <= operand.degree() ? operand[i] : 0.0;
            result[i] = x + y;
        }

        return polynomial { result };
    }

    const polynomial operator* (polynomial operand) const {
        std::vector<double> result (degree() + operand.degree() + 1, 0.0);

        for (int i = 0; i < coefficients.size(); i++) {
            for (int j = 0; j < operand.coefficients.size(); j++) {
                result[i+j] += coefficients[i] * operand.coefficients[j];
            }
        }

        return polynomial { result };
    }

    const polynomial operator- (polynomial operand) const {
        return *this + -operand;
    }

    const polynomial operator- () const {
        return -1.0 * *this;
    }

    friend const polynomial operator* (double x, polynomial operand) {
        for (int i = 0; i < operand.coefficients.size(); i++) {
            operand.coefficients[i] *= x;
        }
        return operand;
    }

    friend const polynomial operator* (polynomial operand, double x) {
        return x * operand;
    }

    friend std::ostream& operator<< (std::ostream& os, const polynomial& p) {
        for (auto c : p.coefficients) {
            os << c << ' ';
        }
        return os;
    }

private:
    std::vector<double> coefficients;
};

std::vector<polynomial> sturm_chain (polynomial p) {
    std::vector<polynomial> chain (p.degree() + 1);

    chain[0] = p;
    chain[1] = p.derivative();

    for (int i = 2; i < chain.size(); i++) {
        chain[i] = -(chain[i-2] % chain[i-1]);
    }

    return chain;
}

#endif
