#pragma once
/**
 * polynomial.h
 * ------------
 * Polynomial  — a callable with explicit monomial coefficients, a std::function
 *               evaluator, and the recurrence metadata that generated it.
 *
 * PolynomialBasis — ordered owning collection of Polynomial objects.
 *
 * Stored field convention
 * -----------------------
 *   p[k].alpha()   = alpha coefficient that built p_k  (alpha_{k-1} in textbooks)
 *   p[k].beta()    = beta  coefficient that built p_k  (beta_{k-1}  in textbooks)
 *   p[k].norm2()   = <p_k, p_k>_w  (true squared L²_w norm)
 *   p[k].scale()   = normalisation factor applied (1.0 if not orthonormal)
 *   p[k].coeffs()  = {c_0, c_1, ..., c_n}  where p_k(x) = c_0 + c_1*x + ... + c_n*x^n
 *
 * Printing
 * --------
 *   p[k].to_string()        — human-readable "a + bx + cx^2 + ..." form
 *   p[k].to_string(digits)  — control number of significant figures
 */

#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace orthopoly {

// ── Single polynomial ─────────────────────────────────────────────────────────

    class Polynomial {
    public:
        /**
         * @param degree   Polynomial degree.
         * @param fn       Callable evaluator.
         * @param coeffs   Monomial coefficients: coeffs[k] is the coefficient of x^k.
         *                 coeffs[0] = constant term, coeffs[degree] = leading coeff.
         * @param alpha    Recurrence alpha stored with this poly.
         * @param beta     Recurrence beta  stored with this poly.
         * @param norm2    Squared L²_w norm.
         * @param scale    Normalisation scale (1.0 if not orthonormal).
         */
        Polynomial(int degree,
                   std::function<double(double)> fn,
                   std::vector<double> coeffs,
                   double alpha = 0.0,
                   double beta  = 0.0,
                   double norm2 = 1.0,
                   double scale = 1.0);

        double operator()(double x) const;
        double eval(double x)       const;

        /// Expose the internal std::function (needed by StieltjesBuilder::gram_matrix).
        const std::function<double(double)>& as_fn_ref() const;

        int                        degree() const;
        double                     alpha()  const;
        double                     beta()   const;
        double                     norm2()  const;
        double                     scale()  const;
        std::string                label()  const;
        const std::vector<double>& coeffs() const;

        /**
         * Format as "a + bx + cx^2 + ..." string.
         *
         * @param digits   Significant figures for each coefficient (default 6).
         * @param eps      Coefficients smaller than eps in absolute value are omitted.
         */
        std::string to_string(int digits = 6, double eps = 1e-12) const;

        /// Print to std::cout (with trailing newline).
        void print(int digits = 6) const;

    private:
        int                           degree_;
        std::function<double(double)> fn_;
        std::vector<double>           coeffs_;   // coeffs_[k] = coeff of x^k
        double                        alpha_, beta_, norm2_, scale_;
    };

// ── Ordered collection ────────────────────────────────────────────────────────

    class PolynomialBasis {
    public:
        PolynomialBasis() = default;

        void push_back(Polynomial p);

        std::size_t size()       const;
        bool        empty()      const;
        int         max_degree() const;

        const Polynomial& operator[](std::size_t i) const;
        Polynomial&       operator[](std::size_t i);

        /// Evaluate all basis polynomials at x.  Returns vector of length size().
        std::vector<double> eval_all(double x) const;

        /// Evaluate all basis polynomials at each x in xs.
        /// Returns matrix[i][j] = p_i(xs[j]).
        std::vector<std::vector<double>> eval_matrix(const std::vector<double>& xs) const;

        /// Print every polynomial in the basis as "a + bx + cx^2 + ...".
        void print_all(int digits = 6) const;

        const std::vector<Polynomial>& polynomials() const;
        std::vector<Polynomial>::const_iterator begin() const;
        std::vector<Polynomial>::const_iterator end()   const;

    private:
        std::vector<Polynomial> polys_;
    };

} // namespace orthopoly