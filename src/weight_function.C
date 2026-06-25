
#include "weight_function.hpp"

#include <cmath>
#include <stdexcept>

namespace orthopoly {

// ── WeightFunction ────────────────────────────────────────────────────────────

    WeightFunction::WeightFunction(double a, double b, std::string name)
            : a_(a), b_(b), name_(std::move(name))
    {
        if (a >= b)
            throw std::invalid_argument("WeightFunction: require a < b");
    }

    std::function<double(double)> WeightFunction::as_fn() const {
        return [this](double x){ return (*this)(x); };
    }

// ── LegendreWeight ────────────────────────────────────────────────────────────

    LegendreWeight::LegendreWeight()
            : WeightFunction(-1.0, 1.0, "Legendre  w(x)=1") {}

    double LegendreWeight::operator()(double /*x*/) const {
        return 1.0;
    }

// ── ChebyshevTWeight ─────────────────────────────────────────────────────────

    ChebyshevTWeight::ChebyshevTWeight()
            : WeightFunction(-1.0, 1.0, "Chebyshev-T  w(x)=1/sqrt(1-x^2)") {}

    double ChebyshevTWeight::operator()(double x) const {
        double v = 1.0 - x * x;
        return (v > 0.0) ? 1.0 / std::sqrt(v) : 0.0;
    }

// ── ChebyshevUWeight ─────────────────────────────────────────────────────────

    ChebyshevUWeight::ChebyshevUWeight()
            : WeightFunction(-1.0, 1.0, "Chebyshev-U  w(x)=sqrt(1-x^2)") {}

    double ChebyshevUWeight::operator()(double x) const {
        double v = 1.0 - x * x;
        return (v > 0.0) ? std::sqrt(v) : 0.0;
    }

// ── HermiteWeight ─────────────────────────────────────────────────────────────

    HermiteWeight::HermiteWeight(double half_range)
            : WeightFunction(-half_range, half_range,
                             "Hermite  w(x)=exp(-x^2)") {}

    double HermiteWeight::operator()(double x) const {
        return std::exp(-x * x);
    }

// ── LaguerreWeight ────────────────────────────────────────────────────────────

    LaguerreWeight::LaguerreWeight(double alpha, double cutoff)
            : WeightFunction(0.0, cutoff,
                             "Laguerre  w(x)=x^" + std::to_string(alpha) + "*exp(-x)"),
              alpha_(alpha)
    {
        if (alpha <= -1.0)
            throw std::invalid_argument("LaguerreWeight: alpha must be > -1");
    }

    double LaguerreWeight::operator()(double x) const {
        if (x <= 0.0) return 0.0;
        return std::pow(x, alpha_) * std::exp(-x);
    }

// ── JacobiWeight ──────────────────────────────────────────────────────────────

    JacobiWeight::JacobiWeight(double alpha, double beta)
            : WeightFunction(-1.0, 1.0,
                             "Jacobi  w(x)=(1-x)^" + std::to_string(alpha)
                             + "*(1+x)^" + std::to_string(beta)),
              alpha_(alpha), beta_(beta)
    {
        if (alpha <= -1.0 || beta <= -1.0)
            throw std::invalid_argument("JacobiWeight: alpha and beta must be > -1");
    }

    double JacobiWeight::operator()(double x) const {
        double lhs = (x <= 1.0)  ? std::pow(1.0 - x, alpha_) : 0.0;
        double rhs = (x >= -1.0) ? std::pow(1.0 + x, beta_)  : 0.0;
        return lhs * rhs;
    }

// ── LambdaWeight ──────────────────────────────────────────────────────────────

    LambdaWeight::LambdaWeight(std::function<double(double)> fn,
                               double a, double b,
                               std::string name)
            : WeightFunction(a, b, std::move(name)), fn_(std::move(fn)) {}

    double LambdaWeight::operator()(double x) const {
        return fn_(x);
    }

// ── Factories ─────────────────────────────────────────────────────────────────

    std::unique_ptr<WeightFunction> make_legendre()
    { return std::make_unique<LegendreWeight>(); }

    std::unique_ptr<WeightFunction> make_chebyshev_T()
    { return std::make_unique<ChebyshevTWeight>(); }

    std::unique_ptr<WeightFunction> make_chebyshev_U()
    { return std::make_unique<ChebyshevUWeight>(); }

    std::unique_ptr<WeightFunction> make_hermite(double half_range)
    { return std::make_unique<HermiteWeight>(half_range); }

    std::unique_ptr<WeightFunction> make_laguerre(double alpha, double cutoff)
    { return std::make_unique<LaguerreWeight>(alpha, cutoff); }

    std::unique_ptr<WeightFunction> make_jacobi(double alpha, double beta)
    { return std::make_unique<JacobiWeight>(alpha, beta); }

} // namespace orthopoly