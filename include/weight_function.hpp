#pragma once
/**
 * weight_function.h
 * -----------------
 * Abstract WeightFunction base class and concrete implementations for all
 * standard 1-D polynomial families.
 *
 * A WeightFunction encapsulates:
 *   - evaluation w(x)
 *   - the integration domain [a, b]
 *   - a human-readable name
 *   LambdaWeight       arbitrary std::function         on any [a, b]
 */

#include <functional>
#include <memory>
#include <string>

namespace orthopoly {

// ── Abstract base ─────────────────────────────────────────────────────────────

    class WeightFunction {
    public:
        WeightFunction(double a, double b, std::string name);
        virtual ~WeightFunction() = default;

        /// Evaluate w(x).
        virtual double operator()(double x) const = 0;

        double             a()    const { return a_; }
        double             b()    const { return b_; }
        const std::string& name() const { return name_; }

        /// Wrap this weight as a std::function (for passing to integration helpers).
        std::function<double(double)> as_fn() const;

    protected:
        double      a_, b_;
        std::string name_;
    };


    class LambdaWeight : public WeightFunction {
    public:
        LambdaWeight(std::function<double(double)> fn,
                     double a, double b,
                     std::string name = "Custom weight");
        double operator()(double x) const override;

    private:
        std::function<double(double)> fn_;
    };

} // namespace orthopoly
