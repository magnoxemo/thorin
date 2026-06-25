#include "orthopoly.hpp"

#include <iomanip>
#include <string>

namespace orthopoly {

    Printer::Printer(std::ostream& out, int col_w)
            : out_(out), col_w_(col_w) {}

    void Printer::banner(const std::string& title) const {
        const int W = 72;
        out_ << "\n" << std::string(W, '=') << "\n"
             << "  " << title << "\n"
             << std::string(W, '=') << "\n";
    }

    void Printer::print_recurrence(const PolynomialBasis& basis) const {
        const int W = 4 + 18 + 18 + 18;
        out_ << "\nRecurrence  p_{n+1} = (x - alpha_n)*p_n - beta_n*p_{n-1}\n";
        out_ << std::string(W, '-') << "\n";
        out_ << std::setw(4)  << "n"
             << std::setw(18) << "alpha_n"
             << std::setw(18) << "beta_n"
             << std::setw(18) << "||p_n||^2_w" << "\n";
        out_ << std::string(W, '-') << "\n";
        out_ << std::fixed << std::setprecision(10);
        for (auto& p : basis)
            out_ << std::setw(4)  << p.degree()
                 << std::setw(18) << p.alpha()
                 << std::setw(18) << p.beta()
                 << std::setw(18) << p.norm2() << "\n";
    }

    void Printer::print_eval_table(const PolynomialBasis& basis,
                                   double a, double b, int n_pts) const
    {
        out_ << "\nPolynomial values  (x from " << a << " to " << b << "):\n";
        out_ << std::setw(10) << "x";
        for (auto& p : basis) out_ << std::setw(col_w_) << p.label();
        out_ << "\n" << std::string(10 + col_w_ * (int)basis.size(), '-') << "\n";

        out_ << std::fixed << std::setprecision(6);
        for (int i = 0; i < n_pts; ++i) {
            double x = (n_pts == 1) ? 0.5 * (a + b)
                                    : a + (b - a) * i / (n_pts - 1);
            out_ << std::setw(10) << x;
            for (auto& p : basis) out_ << std::setw(col_w_) << p(x);
            out_ << "\n";
        }
    }

    void Printer::print_gram_matrix(const std::vector<std::vector<double>>& G,
                                    const std::string& label) const
    {
        const int sz = static_cast<int>(G.size());
        out_ << "\nGram matrix  <p_i, p_j>_w";
        if (!label.empty()) out_ << "  " << label;
        out_ << ":\n";

        out_ << std::setw(6) << "";
        for (int j = 0; j < sz; ++j)
            out_ << std::setw(14) << ("p" + std::to_string(j));
        out_ << "\n";

        out_ << std::scientific << std::setprecision(3);
        for (int i = 0; i < sz; ++i) {
            out_ << std::setw(4) << ("p" + std::to_string(i)) << "  ";
            for (int j = 0; j < sz; ++j) out_ << std::setw(14) << G[i][j];
            out_ << "\n";
        }
    }

} // namespace orthopoly