#include "stieltjes.hpp"

#include <cmath>
#include <stdexcept>

namespace orthopoly {

// ── Constructor / set_quadrature ──────────────────────────────────────────────

    StieltjesBuilder::StieltjesBuilder(const WeightFunction& weight,
                                       BuilderOptions opts)
            : w_(weight), opts_(opts),
              quad_(std::make_unique<GaussLegendreQuadrature>(opts.quad_points))
    {}

    void StieltjesBuilder::set_quadrature(std::unique_ptr<QuadratureRule> q) {
        quad_ = std::move(q);
    }

    double StieltjesBuilder::ip(const ScalarFn& f, const ScalarFn& g,
                                const ScalarFn& wfn, double a, double b) const
    {
        return inner_product(f, g, wfn, a, b, *quad_);
    }

// ── Coefficient arithmetic ────────────────────────────────────────────────────
// All operations on coefficient vectors c where c[k] = coefficient of x^k.

// Multiply coefficient vector by scalar
    static std::vector<double> scale_coeffs(const std::vector<double>& c, double s) {
        std::vector<double> r(c.size());
        for (std::size_t i = 0; i < c.size(); ++i) r[i] = c[i] * s;
        return r;
    }

// c1 - c2  (vectors may differ in length)
    static std::vector<double> sub_coeffs(const std::vector<double>& c1,
                                          const std::vector<double>& c2)
    {
        std::size_t n = std::max(c1.size(), c2.size());
        std::vector<double> r(n, 0.0);
        for (std::size_t i = 0; i < c1.size(); ++i) r[i] += c1[i];
        for (std::size_t i = 0; i < c2.size(); ++i) r[i] -= c2[i];
        return r;
    }

// Multiply coefficient vector by x  (shift all coefficients up by one degree)
    static std::vector<double> mul_x(const std::vector<double>& c) {
        std::vector<double> r(c.size() + 1, 0.0);
        for (std::size_t i = 0; i < c.size(); ++i) r[i + 1] = c[i];
        return r;
    }

// Multiply coefficient vector by (-alpha) and add to x*c
// i.e. compute (x - alpha) * c
    static std::vector<double> mul_x_minus_alpha(const std::vector<double>& c,
                                                 double alpha)
    {
        // (x - alpha)*p = x*p - alpha*p
        auto xp = mul_x(c);
        for (std::size_t i = 0; i < c.size(); ++i) xp[i] -= alpha * c[i];
        return xp;
    }

// ── build ─────────────────────────────────────────────────────────────────────

    PolynomialBasis StieltjesBuilder::build(int N) const {
        if (N < 0)
            throw std::invalid_argument("StieltjesBuilder::build: N must be >= 0");

        const double   a   = w_.a();
        const double   b   = w_.b();
        const ScalarFn wfn = w_.as_fn();

        PolynomialBasis basis;

        // ── p_0 = 1 ──────────────────────────────────────────────────────────────
        ScalarFn            p0_raw    = [](double) { return 1.0; };
        std::vector<double> p0_coeffs = {1.0};

        double n0 = ip(p0_raw, p0_raw, wfn, a, b);
        if (!std::isfinite(n0) || n0 < opts_.zero_norm_tol)
            throw std::runtime_error(
                    "StieltjesBuilder: integral of w(x) is not positive — "
                    "check weight function and interval.");

        double   scale0  = opts_.orthonormal ? 1.0 / std::sqrt(n0) : 1.0;
        ScalarFn stored0 = opts_.orthonormal
                           ? ScalarFn([scale0](double) { return scale0; })
                           : p0_raw;
        auto coeffs0 = opts_.orthonormal ? scale_coeffs(p0_coeffs, scale0) : p0_coeffs;

        basis.push_back(Polynomial(0, stored0, coeffs0, 0.0, 0.0, n0, scale0));

        if (N == 0) return basis;

        // ── Recurrence for n = 0, 1, ..., N-1 ────────────────────────────────────
        std::vector<ScalarFn>            raw;
        std::vector<double>              raw_norm2;
        std::vector<std::vector<double>> raw_coeffs;  // exact coefficients

        raw.push_back(p0_raw);
        raw_norm2.push_back(n0);
        raw_coeffs.push_back(p0_coeffs);

        for (int n = 0; n < N; ++n) {
            const ScalarFn&            pn   = raw[n];
            double                     nn   = raw_norm2[n];
            const std::vector<double>& cn   = raw_coeffs[n];

            // alpha_n = <x*p_n, p_n>_w / <p_n, p_n>_w
            ScalarFn xpn  = [&pn](double x){ return x * pn(x); };
            double   alpha = ip(xpn, pn, wfn, a, b) / nn;

            // beta_n = <p_n, p_n>_w / <p_{n-1}, p_{n-1}>_w   (0 for n=0)
            double beta = (n == 0) ? 0.0 : nn / raw_norm2[n - 1];

            // Evaluator: p_{n+1}(x) = (x - alpha)*p_n(x) - beta*p_{n-1}(x)
            ScalarFn pn_cap   = pn;
            ScalarFn pnm1_cap = (n == 0)
                                ? ScalarFn([](double) { return 0.0; })
                                : raw[n - 1];
            double alpha_c = alpha, beta_c = beta;

            ScalarFn pn1_raw = [pn_cap, pnm1_cap, alpha_c, beta_c](double x) {
                return (x - alpha_c) * pn_cap(x) - beta_c * pnm1_cap(x);
            };

            // Exact coefficients: c_{n+1} = (x - alpha)*c_n - beta*c_{n-1}
            auto cn1 = mul_x_minus_alpha(cn, alpha);
            if (n > 0) {
                auto beta_cn_1 = scale_coeffs(raw_coeffs[n - 1], beta);
                cn1 = sub_coeffs(cn1, beta_cn_1);
            }

            double nn1 = ip(pn1_raw, pn1_raw, wfn, a, b);
            if (!std::isfinite(nn1) || nn1 < opts_.zero_norm_tol)
                throw std::runtime_error(
                        "StieltjesBuilder: near-zero norm at degree "
                        + std::to_string(n + 1)
                        + " — try more quadrature points.");

            raw.push_back(pn1_raw);
            raw_norm2.push_back(nn1);
            raw_coeffs.push_back(cn1);

            // Build stored evaluator (optionally normalised)
            double   scaleN = opts_.orthonormal ? 1.0 / std::sqrt(nn1) : 1.0;
            ScalarFn stored = opts_.orthonormal
                              ? ScalarFn([pn1_raw, scaleN](double x){ return scaleN * pn1_raw(x); })
                              : pn1_raw;
            auto stored_coeffs = opts_.orthonormal ? scale_coeffs(cn1, scaleN) : cn1;

            basis.push_back(Polynomial(n + 1, stored, stored_coeffs,
                                       alpha, beta, nn1, scaleN));
        }

        return basis;
    }

// ── gram_matrix ───────────────────────────────────────────────────────────────

    std::vector<std::vector<double>>
    StieltjesBuilder::gram_matrix(const PolynomialBasis& basis) const {
        const int    sz  = static_cast<int>(basis.size());
        const double a   = w_.a();
        const double b   = w_.b();
        const ScalarFn wfn = w_.as_fn();

        std::vector<std::vector<double>> G(sz, std::vector<double>(sz, 0.0));
        for (int i = 0; i < sz; ++i)
            for (int j = i; j < sz; ++j) {
                double v = ip(basis[i].as_fn_ref(), basis[j].as_fn_ref(), wfn, a, b);
                G[i][j] = v;
                G[j][i] = v;
            }
        return G;
    }

} // namespace orthopoly