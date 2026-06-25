#pragma once
/**
 * orthopoly.h
 * -----------
 * Single header that pulls in the entire orthopoly framework.
 * Also provides a Printer utility for formatted console output.
 *
 * Usage:
 *   #include "orthopoly.h"
 *   using namespace orthopoly;
 */

#include "polynomial.hpp"
#include "quadrature.hpp"
#include "stieltjes.hpp"
#include "weight_function.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace orthopoly {

// ── Pretty printer ────────────────────────────────────────────────────────────

    class Printer {
    public:
        explicit Printer(std::ostream& out = std::cout, int col_w = 16);

        /// Print a section banner.
        void banner(const std::string& title) const;

        /// Print recurrence coefficients alpha_n, beta_n, ||p_n||^2 for each polynomial.
        void print_recurrence(const PolynomialBasis& basis) const;

        /// Print a table of polynomial values at n_pts evenly-spaced points in [a, b].
        void print_eval_table(const PolynomialBasis& basis,
                              double a, double b,
                              int n_pts = 7) const;

        /// Print the full Gram matrix.
        void print_gram_matrix(const std::vector<std::vector<double>>& G,
                               const std::string& label = "") const;

    private:
        std::ostream& out_;
        int           col_w_;
    };

} // namespace orthopoly