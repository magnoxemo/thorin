#include "orthopoly.hpp"

using namespace orthopoly;

int main() {
    std::function<double(double )> f ;
    LambdaWeight w([](double x){ return (1.0+x)*(1-x)+(1-x); }, -1.0, 1.0, "Cauchy");
    auto builder = StieltjesBuilder(w);
    builder.opts_.quad_points = 10000;
    auto basis = builder.build(10);
    basis.print_all();
    Printer p;
    p.print_gram_matrix(builder.gram_matrix(basis));
    return 0;
}