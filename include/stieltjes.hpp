#pragma once
/**
 * stieltjes.h
 * -----------
 * StieltjesBuilder: generates orthogonal or orthonormal polynomial bases
 * for an arbitrary weight function using the three-term recurrence:
 *
 *   p_{n+1}(x) = (x - alpha_n) * p_n(x) - beta_n * p_{n-1}(x)
 *
 * where
 *   alpha_n = <x p_n, p_n>_w / <p_n, p_n>_w
 *   beta_n  = <p_n, p_n>_w   / <p_{n-1}, p_{n-1}>_w   (beta_0 = 0)
 *
 * Inner products are computed via a QuadratureRule (default: 80-point GL).
 * For nearly-singular weights, increase quad_points in BuilderOptions.
 */

#include "polynomial.hpp"
#include "quadrature.hpp"
#include "weight_function.hpp"

#include <memory>

namespace orthopoly {

// ── Options ───────────────────────────────────────────────────────────────────

    struct BuilderOptions {
        bool   orthonormal    = false;  ///< normalise so <p_i, p_i>_w = 1
        int    quad_points    = 80;     ///< Gauss-Legendre nodes for inner products
        double zero_norm_tol  = 1e-14; ///< abort if a squared norm falls below this
    };

// ── Builder ───────────────────────────────────────────────────────────────────

    class StieltjesBuilder {
    public:
        /**
         * @param weight  Weight function (borrowed; must outlive the builder).
         * @param opts    Tuning options.
         */
        explicit StieltjesBuilder(const WeightFunction& weight,
                                  BuilderOptions opts = {});

        /// Replace the quadrature rule (e.g. for specialised nodes near singularities).
        void set_quadrature(std::unique_ptr<QuadratureRule> q);

        /**
         * Build and return the polynomial basis p_0, ..., p_N.
         * Throws std::runtime_error if any squared norm is non-positive.
         */
        PolynomialBasis build(int N) const;

        /**
         * Compute the Gram matrix G[i][j] = <p_i, p_j>_w for a given basis.
         * For an orthogonal basis this is diagonal; for orthonormal it is the identity.
         */
        std::vector<std::vector<double>> gram_matrix(const PolynomialBasis& basis) const;

        BuilderOptions                  opts_;
    private:
        double ip(const ScalarFn& f, const ScalarFn& g,
                  const ScalarFn& wfn, double a, double b) const;

        const WeightFunction&           w_;
        std::unique_ptr<QuadratureRule> quad_;
    };

} // namespace orthopoly