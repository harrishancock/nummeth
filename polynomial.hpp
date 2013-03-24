#ifndef NUMMETH_POLYNOMIAL_HPP
#define NUMMETH_POLYNOMIAL_HPP

#include <vector>
#include <iostream>

class polynomial {
public:
    explicit polynomial (std::vector<double> cs) : coefficients(cs) { }

    double operator() (const double x) const {
        double p = 0;
        for (auto c : coefficients) {
            p = p * x + c;
        }
        return p;
    }

    double operator[] (const ssize_t i) const {
        return coefficients[i];
    }

    size_t size () const {
        return coefficients.size();
    }

    polynomial derivative () const {
        if (coefficients.size() <= 1) {
            return polynomial { std::vector<double>() };
        }

        std::vector<double> cs { coefficients.begin(), coefficients.end() - 1 };
        for (int i = 0; i < cs.size(); i++) {
            cs[i] *= cs.size() - i;
        }
        
        return polynomial { cs };
    }

    std::pair<polynomial, polynomial> operator/ (polynomial divisor) const {
        int rsize = divisor.size() - 1;
        int qsize = coefficients.size() - rsize;

        std::vector<double> remainder { coefficients.begin() + qsize, coefficients.end() };
        std::vector<double> quotient { coefficients.begin(), coefficients.begin() + qsize };

        /* FIXME rewrite this more functionally, with lists and shit */

        for (int i = 0; i < quotient.size(); i++) {
            for (int j = 1; j < divisor.size(); j++) {
                if (i - j >= 0) {
                    quotient[i] -= divisor[j] * quotient[i-j];
                }
            }
            quotient[i] /= divisor[0];
        }

        for (int i = 0; i < remainder.size(); i++) {
            for (int j = i+1; j < divisor.size(); j++) {
                if (qsize + i - j >= 0) {
                    remainder[i] -= divisor[j] * quotient[qsize + i - j];
                }
            }
        }

        return { polynomial { quotient }, polynomial { remainder } };
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

#endif
