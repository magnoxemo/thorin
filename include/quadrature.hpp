#pragma once
/**
 * quadrature.h
 * ------------
 * Abstract numerical integration rule and Gauss-Legendre implementation.
 *
 * QuadratureRule  — base class; owns pre-computed nodes/weights on a reference
 *                   domain; integrates arbitrary functions over [a, b].
 * GaussLegendreQuadrature — n-point GL rule, exact for polynomials of degree ≤ 2n-1.
 */

#include <functional>
#include <vector>

namespace orthopoly {

    using ScalarFn = std::function<double(double)>;

// ── Abstract base ─────────────────────────────────────────────────────────────

    class QuadratureRule {
    public:
        explicit QuadratureRule(int n_points);
        virtual ~QuadratureRule() = default;

        /// Integrate f over [a, b] using the stored nodes/weights.
        double integrate(const ScalarFn& f, double a, double b) const;

        int    n_points()              const { return n_; }
        double node  (int k)           const { return xi_[k]; }
        double weight(int k)           const { return w_[k];  }
        const std::vector<double>& nodes()   const { return xi_; }
        const std::vector<double>& weights() const { return w_;  }

    protected:
        int                 n_;
        std::vector<double> xi_;   // reference-domain nodes
        std::vector<double> w_;    // quadrature weights
    };

// ── Gauss-Legendre on [-1, 1] ─────────────────────────────────────────────────

    class GaussLegendreQuadrature : public QuadratureRule {
    public:
        explicit GaussLegendreQuadrature(int n);

    private:
        void compute_nodes_weights();
    };

// ── Free helper ───────────────────────────────────────────────────────────────

/// Compute <f, g>_w = ∫_a^b w(x) f(x) g(x) dx via the given quadrature rule.
    double inner_product(const ScalarFn& f,
                         const ScalarFn& g,
                         const ScalarFn& w,
                         double a, double b,
                         const QuadratureRule& quad);

} // namespace orthopoly