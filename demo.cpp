#include "orthopoly.hpp"

#include <cmath>

using namespace orthopoly;

int main() {

    LambdaWeight w([](double x){ return std::abs(1.0+x-3*x*x+sin(x)); }, -1.0, 1.0, "random");
    auto builder = StieltjesBuilder(w);
    builder.options().quad_points = 10000;
    auto basis = builder.build(20);
    basis.print_all();

    Printer().print_gram_matrix(builder.gram_matrix(basis));
    return 0;
}