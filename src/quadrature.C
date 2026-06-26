#include "quadrature.hpp"

#include <cmath>
#include <stdexcept>

namespace orthopoly {

// ── QuadratureRule ────────────────────────────────────────────────────────────

    QuadratureRule::QuadratureRule(int n_points) : n_(n_points) {
        if (n_points < 1)
            throw std::invalid_argument("QuadratureRule: n_points must be >= 1");
    }

    double QuadratureRule::integrate(const ScalarFn& f, double a, double b) const {
        if (a >= b)
            throw std::invalid_argument("QuadratureRule::integrate: require a < b");
        double mid  = 0.5 * (a + b);
        double half = 0.5 * (b - a);
        double sum  = 0.0;
        for (int k = 0; k < n_; ++k)
            sum += w_[k] * f(mid + half * xi_[k]);
        return sum * half;
    }

// ── GaussLegendreQuadrature ───────────────────────────────────────────────────

    GaussLegendreQuadrature::GaussLegendreQuadrature(int n)
            : QuadratureRule(n)
    {
        xi_.resize(n);
        w_.resize(n);
        compute_nodes_weights();
    }

    void GaussLegendreQuadrature::compute_nodes_weights() {
        const int half = (n_ + 1) / 2;
        for (int i = 0; i < half; ++i) {
            // Asymptotic initial guess for the i-th root of P_n
            double xi = std::cos(M_PI * (i + 0.75) / (n_ + 0.5));
            double p1, p2, p3, dp;
            // Newton refinement on P_n
            for (int iter = 0; iter < 1000; ++iter)
            {
                p1 = 1.0; p2 = 0.0;

                for (int j = 1; j <= n_; ++j)
                {
                    p3 = p2; p2 = p1;
                    p1 = ((2*j - 1) * xi * p2 - (j - 1) * p3) / j;
                }

                // Derivative via: (1-x²) P_n'(x) = n (P_{n-1}(x) - x P_n(x))
                dp = n_ * (xi * p1 - p2) / (xi * xi - 1.0);
                double dxi = p1 / dp;
                xi -= dxi;
                if (std::abs(dxi) < 1e-20) break;
            }
            // Symmetric placement
            xi_[i]          = -xi;
            xi_[n_ - 1 - i] =  xi;
            double w_val = 2.0 / ((1.0 - xi * xi) * dp * dp);
            w_[i]          = w_val;
            w_[n_ - 1 - i] = w_val;
        }
    }

// ── Free helper ───────────────────────────────────────────────────────────────

    double inner_product(const ScalarFn& f,
                         const ScalarFn& g,
                         const ScalarFn& w,
                         double a, double b,
                         const QuadratureRule& quad)
    {
        return quad.integrate([&](double x){ return w(x) * f(x) * g(x); }, a, b);
    }

} // namespace orthopoly