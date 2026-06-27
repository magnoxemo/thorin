
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


// ── LambdaWeight ──────────────────────────────────────────────────────────────

    LambdaWeight::LambdaWeight(std::function<double(double)> fn,
                               double a, double b,
                               std::string name)
            : WeightFunction(a, b, std::move(name)), fn_(std::move(fn)) {}

    double LambdaWeight::operator()(double x) const {
        return fn_(x);
    }


} // namespace orthopoly