#include "weight_function.hpp"


namespace tests{
    
    class LegendreWeight : public orthopoly::WeightFunction {
    public:
        LegendreWeight();
        double operator()(double x) const override;
    };

    class ChebyshevTWeight : public orthopoly::WeightFunction {
    public:
        ChebyshevTWeight();
        double operator()(double x) const override;
    };

    class ChebyshevUWeight : public orthopoly::WeightFunction {
    public:
        ChebyshevUWeight();
        double operator()(double x) const override;
    };

    class HermiteWeight : public orthopoly::WeightFunction {
    public:
        /// @param half_range  domain is [-half_range, +half_range]
        explicit HermiteWeight(double half_range = 6.0);
        double operator()(double x) const override;
    };

    class LaguerreWeight : public orthopoly::WeightFunction {
    public:
        /// @param alpha   generalized parameter (must be > -1)
        /// @param cutoff  upper integration limit
        explicit LaguerreWeight(double alpha = 0.0, double cutoff = 30.0);
        double operator()(double x) const override;
        double alpha() const { return alpha_; }

    private:
        double alpha_;
    };

    class JacobiWeight : public orthopoly::WeightFunction {
    public:
        /// @param alpha, beta  must be > -1
        JacobiWeight(double alpha, double beta);
        double operator()(double x) const override;
        double alpha() const { return alpha_; }
        double beta()  const { return beta_;  }

    private:
        double alpha_, beta_;
    };


// ── Factory helpers ───────────────────────────────────────────────────────────

    std::unique_ptr<orthopoly::WeightFunction> make_legendre();
    std::unique_ptr<orthopoly::WeightFunction> make_chebyshev_T();
    std::unique_ptr<orthopoly::WeightFunction> make_chebyshev_U();
    std::unique_ptr<orthopoly::WeightFunction> make_hermite(double half_range = 6.0);
    std::unique_ptr<orthopoly::WeightFunction> make_laguerre(double alpha = 0.0, double cutoff = 30.0);
    std::unique_ptr<orthopoly::WeightFunction> make_jacobi(double alpha, double beta);

}