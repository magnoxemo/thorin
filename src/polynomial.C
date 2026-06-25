#include "polynomial.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace orthopoly {

// ── Polynomial ────────────────────────────────────────────────────────────────

    Polynomial::Polynomial(int degree,
                           std::function<double(double)> fn,
                           std::vector<double> coeffs,
                           double alpha,
                           double beta,
                           double norm2,
                           double scale)
            : degree_(degree), fn_(std::move(fn)), coeffs_(std::move(coeffs)),
              alpha_(alpha), beta_(beta), norm2_(norm2), scale_(scale)
    {}

    double Polynomial::operator()(double x) const { return fn_(x); }
    double Polynomial::eval(double x)       const { return fn_(x); }

    const std::function<double(double)>& Polynomial::as_fn_ref() const { return fn_; }

    int                        Polynomial::degree() const { return degree_; }
    double                     Polynomial::alpha()  const { return alpha_;  }
    double                     Polynomial::beta()   const { return beta_;   }
    double                     Polynomial::norm2()  const { return norm2_;  }
    double                     Polynomial::scale()  const { return scale_;  }
    const std::vector<double>& Polynomial::coeffs() const { return coeffs_; }
    std::string                Polynomial::label()  const {
        return "p" + std::to_string(degree_);
    }

// ── to_string ─────────────────────────────────────────────────────────────────

    std::string Polynomial::to_string(int digits, double eps) const {
        // Format one coefficient value with the requested significant figures
        auto fmt = [&](double c) -> std::string {
            std::ostringstream ss;
            ss << std::setprecision(digits) << c;
            return ss.str();
        };

        std::string result;
        bool first = true;

        for (int k = 0; k <= degree_; ++k) {
            double c = coeffs_[k];
            if (std::abs(c) < eps) continue;   // skip near-zero terms

            // Sign / separator
            if (first) {
                if (c < 0) result += "-";
            } else {
                result += (c < 0 ? " - " : " + ");
            }
            double ac = std::abs(c);

            // Coefficient — omit "1" when multiplied by a power of x (except x^0)
            bool is_one = (std::abs(ac - 1.0) < eps);

            if (k == 0) {
                // Constant term: always print the number
                result += fmt(ac);
            } else if (k == 1) {
                // Linear term
                if (!is_one) result += fmt(ac);
                result += "x";
            } else {
                // Higher powers
                if (!is_one) result += fmt(ac);
                result += "x^" + std::to_string(k);
            }

            first = false;
        }

        if (result.empty()) result = "0";
        return result;
    }

    void Polynomial::print(int digits) const {
        std::cout << label() << "(x) = " << to_string(digits) << "\n";
    }

// ── PolynomialBasis ───────────────────────────────────────────────────────────

    void PolynomialBasis::push_back(Polynomial p) { polys_.push_back(std::move(p)); }

    std::size_t PolynomialBasis::size()  const { return polys_.size();  }
    bool        PolynomialBasis::empty() const { return polys_.empty(); }

    int PolynomialBasis::max_degree() const {
        return polys_.empty() ? -1 : polys_.back().degree();
    }

    const Polynomial& PolynomialBasis::operator[](std::size_t i) const {
        if (i >= polys_.size())
            throw std::out_of_range("PolynomialBasis: index out of range");
        return polys_[i];
    }

    Polynomial& PolynomialBasis::operator[](std::size_t i) {
        if (i >= polys_.size())
            throw std::out_of_range("PolynomialBasis: index out of range");
        return polys_[i];
    }

    std::vector<double> PolynomialBasis::eval_all(double x) const {
        std::vector<double> vals;
        vals.reserve(polys_.size());
        for (auto& p : polys_) vals.push_back(p(x));
        return vals;
    }

    std::vector<std::vector<double>>
    PolynomialBasis::eval_matrix(const std::vector<double>& xs) const {
        std::vector<std::vector<double>> mat(polys_.size(),
                                             std::vector<double>(xs.size()));
        for (std::size_t i = 0; i < polys_.size(); ++i)
            for (std::size_t j = 0; j < xs.size(); ++j)
                mat[i][j] = polys_[i](xs[j]);
        return mat;
    }

    void PolynomialBasis::print_all(int digits) const {
        for (auto& p : polys_) p.print(digits);
    }

    const std::vector<Polynomial>& PolynomialBasis::polynomials() const { return polys_; }

    std::vector<Polynomial>::const_iterator PolynomialBasis::begin() const { return polys_.begin(); }
    std::vector<Polynomial>::const_iterator PolynomialBasis::end()   const { return polys_.end();   }

} // namespace orthopoly