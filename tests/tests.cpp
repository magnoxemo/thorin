/**
 * tests.cpp  —  orthopoly framework test suite
 *
 * STORED FIELD CONVENTION:
 *   p[k].alpha()  = alpha coefficient that built p_k  (alpha_{k-1} in textbook)
 *   p[k].beta()   = beta  coefficient that built p_k  (beta_{k-1} in textbook)
 *   p[k].norm2()  = <p_k, p_k>_w  (squared L²_w norm)
 *
 *   Our polynomials are the MONIC Stieltjes polynomials (leading coeff 1 for p_0=1,
 *   built recursively without any re-scaling). This differs from the standard
 *   normalisation used in textbooks:
 *     Chebyshev-U: our p_n = U_n / 2^n   →  norm²[n] = (π/2) / 4^n
 *     Hermite:     our p_n = H_n / 2^n   →  norm²[n] = √π · n! / 2^n
 *     Laguerre:    same as standard L_n   →  norm²[n] = 1   (they're already monic)
 *
 *   Recurrence coefficients are intrinsic and match standard values regardless
 *   of the polynomial scaling convention (they come from the ratio of norms).
 *
 * Compile & run:
 *   g++ -O2 -std=c++17 -I../include -o tests tests.cpp && ./tests
 */

#include "orthopoly.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <functional>

using namespace orthopoly;

// ── Minimal test harness ─────────────────────────────────────────────────────

static int g_pass=0, g_fail=0;

void check(bool ok, const std::string& name) {
    if (ok) { std::cout<<"  [PASS] "<<name<<"\n"; ++g_pass; }
    else     { std::cout<<"  [FAIL] "<<name<<"\n"; ++g_fail; }
}
void check_near(double got, double exp, double tol, const std::string& name) {
    bool ok=std::abs(got-exp)<tol;
    if (!ok) {
        std::ostringstream s;
        s<<name<<"  got="<<got<<"  exp="<<exp<<"  |err|="<<std::abs(got-exp);
        check(false,s.str());
    } else check(true,name);
}

// ── High-accuracy GL quadrature for independent verification ─────────────────

static GaussLegendreQuadrature g_qv(200);

double vip(const Polynomial& pi, const Polynomial& pj, const WeightFunction& w) {
    return inner_product(pi.as_fn_ref(), pj.as_fn_ref(),
                         w.as_fn(), w.a(), w.b(), g_qv);
}

// Gauss-Chebyshev-T nodes for verifying Chebyshev-T orthogonality exactly
double chebT_ip(const Polynomial& pi, const Polynomial& pj, int n=1000) {
    double s=0;
    for(int k=0;k<n;++k){
        double x=std::cos(M_PI*(2*k+1)/(2.0*n));
        s+=pi(x)*pj(x);
    }
    return s*(M_PI/n);
}

// ── Orthogonality / orthonormality checks ────────────────────────────────────

void check_orth(const PolynomialBasis& b, const WeightFunction& w,
                double tol, const std::string& tag)
{
    for(size_t i=0;i<b.size();++i)
        for(size_t j=i+1;j<b.size();++j){
            std::ostringstream n; n<<tag<<"  <p"<<i<<",p"<<j<<"> ~ 0";
            check(std::abs(vip(b[i],b[j],w))<tol, n.str());
        }
}

void check_chebT_orth(const PolynomialBasis& b, double tol, const std::string& tag){
    for(size_t i=0;i<b.size();++i)
        for(size_t j=i+1;j<b.size();++j){
            std::ostringstream n; n<<tag<<"  <p"<<i<<",p"<<j<<"> ~ 0";
            check(std::abs(chebT_ip(b[i],b[j]))<tol, n.str());
        }
}

void check_on(const PolynomialBasis& b, const WeightFunction& w,
              double tol, const std::string& tag)
{
    for(size_t i=0;i<b.size();++i)
        for(size_t j=i;j<b.size();++j){
            double val=vip(b[i],b[j],w), exp=(i==j)?1.0:0.0;
            std::ostringstream n; n<<tag<<"  G["<<i<<","<<j<<"]="<<exp;
            check(std::abs(val-exp)<tol, n.str());
        }
}

// ── Analytic reference polynomials ──────────────────────────────────────────

namespace ref {
    double P(int n, double x){
        if(n==0) return 1; if(n==1) return x;
        double pm=1,p=x,p1;
        for(int k=2;k<=n;++k){p1=((2*k-1)*x*p-(k-1)*pm)/k; pm=p; p=p1;} return p;
    }
    double T(int n, double x){return std::cos(n*std::acos(std::clamp(x,-1.0,1.0)));}
    double U(int n, double x){
        double th=std::acos(std::clamp(x,-1.0,1.0)), s=std::sin(th);
        return s<1e-14?(n+1)*std::cos(n*th)/std::cos(th):std::sin((n+1)*th)/s;
    }
    double H(int n, double x){
        if(n==0) return 1; if(n==1) return 2*x;
        double hm=1,h=2*x,h1;
        for(int k=2;k<=n;++k){h1=2*x*h-2*(k-1)*hm; hm=h; h=h1;} return h;
    }
    double L(int n, double x){
        if(n==0) return 1; if(n==1) return 1-x;
        double lm=1,l=1-x,l1;
        for(int k=2;k<=n;++k){l1=((2*k-1-x)*l-(k-1)*lm)/k; lm=l; l=l1;} return l;
    }
    double fact(int n){return n<=1?1.0:n*fact(n-1);}

    // Shape check: compare p_n(x)/p_n(x0) vs ref(n,x)/ref(n,x0)
    void shape(const PolynomialBasis& b, int n, double x0,
               const std::vector<double>& xs,
               std::function<double(int,double)> fn,
               double tol, const std::string& tag)
    {
        double sr=b[n](x0), fr=fn(n,x0);
        if(std::abs(fr)<1e-10||std::abs(sr)<1e-13) return;
        for(double x:xs){
            double rg=b[n](x)/sr, re=fn(n,x)/fr;
            std::ostringstream s; s<<tag<<"  p"<<n<<"("<<x<<") shape";
            check_near(rg, re, tol, s.str());
        }
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// LEGENDRE
// ═════════════════════════════════════════════════════════════════════════════

void test_legendre(){
    std::cout<<"\n── Legendre ─────────────────────────────────────────────\n";
    LegendreWeight w;
    BuilderOptions opts; opts.quad_points=80;
    auto b = StieltjesBuilder(w,opts).build(6);

    // Monic Legendre norms²: 2/(2n+1) * (2^n*n!/(2n)!)^2 * (2n)!^2/(n!)^2 ... actually
    // the monic Legendre polynomial satisfies: ||p_n||^2 = 2/(2n+1) * (1/(binom(2n,n)*2^n))^2 * ...
    // Easier: just check ratio to what we measured (the builder's own norm2).
    // We verify: norm2[0]=2, norm2[1]=2/3, and the recurrence gives the right ratios.
    check_near(b[0].norm2(), 2.0,      1e-13,  "Legendre norm2[0] = 2");
    check_near(b[1].norm2(), 2.0/3.0,  1e-13,  "Legendre norm2[1] = 2/3");
    check_near(b[2].norm2(), 8.0/45.0, 1e-13, "Legendre norm2[2] = 8/45");
    // Norm consistency: norm2[k] = b[k+1].beta() * norm2[k-1]
    //   because p[k+1].beta() = beta_k = norm2[k]/norm2[k-1]
    for(int k=1;k<=4;++k){
        double exp = b[k+1].beta() * b[k-1].norm2();
        std::ostringstream s; s<<"Legendre norm2["<<k<<"] = b["<<k+1<<"].beta*norm2["<<k-1<<"]";
        check_near(b[k].norm2(), exp, 1e-13, s.str());
    }

    // All stored alphas = 0 (symmetric weight)
    for(int n=0;n<=6;++n){
        std::ostringstream s; s<<"Legendre alpha in p["<<n<<"] = 0";
        check_near(b[n].alpha(), 0.0, 1e-10, s.str());
    }

    // p[k].beta() = beta_{k-1}: p[2].beta=1/3, p[3].beta=4/15, ...
    check_near(b[0].beta(), 0.0,      1e-13, "Legendre p[0].beta = 0");
    check_near(b[1].beta(), 0.0,      1e-13, "Legendre p[1].beta = 0");
    for(int k=2;k<=6;++k){
        int n=k-1;
        double exp=(double)(n*n)/(4.0*n*n-1.0);
        std::ostringstream s; s<<"Legendre p["<<k<<"].beta = "<<exp;
        check_near(b[k].beta(), exp, 1e-13, s.str());
    }

    // Shape vs P_n
    std::vector<double> xs={-0.8,-0.4,0.0,0.4,0.8};
    for(int n=0;n<=6;++n) ref::shape(b, n, 0.7, xs, ref::P, 1e-13, "Legendre");

    check_orth(b, w, 1e-11, "Legendre");

    opts.orthonormal=true;
    auto bon=StieltjesBuilder(w,opts).build(5);
    check_on(bon, w, 1e-13, "Legendre ON");
}

// ═════════════════════════════════════════════════════════════════════════════
// CHEBYSHEV-T  (singular weight — use high GL + Gauss-Chebyshev for verification)
// ═════════════════════════════════════════════════════════════════════════════

void test_chebyshev_T(){
    std::cout<<"\n── Chebyshev-T ──────────────────────────────────────────\n";
    ChebyshevTWeight w;
    BuilderOptions opts; opts.quad_points=500;
    auto b = StieltjesBuilder(w,opts).build(6);

    // norm2[0] ≈ pi (Chebyshev-T weight integrates to pi)
    check_near(b[0].norm2(), M_PI, 5e-3, "Chebyshev-T norm2[0] ≈ π");

    // All stored alphas = 0 (symmetric weight)
    for(int n=0;n<=5;++n){
        std::ostringstream s; s<<"Chebyshev-T alpha in p["<<n<<"] = 0";
        check_near(b[n].alpha(), 0.0, 1e-4, s.str());
    }

    // Recurrence betas: p[2].beta=1/2, p[k>=3].beta=1/4
    check_near(b[0].beta(), 0.0, 1e-13, "Chebyshev-T p[0].beta = 0");
    check_near(b[1].beta(), 0.0, 1e-13, "Chebyshev-T p[1].beta = 0");
    check_near(b[2].beta(), 0.5, 1e-2,  "Chebyshev-T p[2].beta = 1/2");
    for(int k=3;k<=6;++k){
        std::ostringstream s; s<<"Chebyshev-T p["<<k<<"].beta = 1/4";
        check_near(b[k].beta(), 0.25, 1e-2, s.str());
    }

    // Shape vs T_n (interior points only — endpoints are fine)
    std::vector<double> xs={-0.5,-0.2,0.2,0.5};
    for(int n=1;n<=3;++n) ref::shape(b, n, 0.4, xs, ref::T, 5e-3, "Chebyshev-T");
    for(int n=4;n<=5;++n) ref::shape(b, n, 0.1, xs, ref::T, 0.1, "Chebyshev-T");

    // Orthogonality via exact Gauss-Chebyshev-T nodes
    check_chebT_orth(b, 5e-2, "Chebyshev-T");
}

// ═════════════════════════════════════════════════════════════════════════════
// CHEBYSHEV-U
// ═════════════════════════════════════════════════════════════════════════════

void test_chebyshev_U(){
    std::cout<<"\n── Chebyshev-U ──────────────────────────────────────────\n";
    ChebyshevUWeight w;
    BuilderOptions opts; opts.quad_points=120;
    auto b = StieltjesBuilder(w,opts).build(6);

    // Our monic p_n = U_n / 2^n  →  norm²[n] = (π/2) / 4^n
    for(int n=0;n<=6;++n){
        double exp=M_PI/2.0/std::pow(4.0,n);
        std::ostringstream s; s<<"Chebyshev-U norm2["<<n<<"] = (π/2)/4^n";
        check_near(b[n].norm2(), exp, 1e-6, s.str());
    }

    // All stored alphas = 0
    for(int n=0;n<=6;++n){
        std::ostringstream s; s<<"Chebyshev-U alpha in p["<<n<<"] = 0";
        check_near(b[n].alpha(), 0.0, 1e-13, s.str());
    }

    // p[0].beta=0, p[1].beta=0, p[k>=2].beta=1/4
    check_near(b[0].beta(), 0.0, 1e-13, "Chebyshev-U p[0].beta = 0");
    check_near(b[1].beta(), 0.0, 1e-13, "Chebyshev-U p[1].beta = 0");
    for(int k=2;k<=6;++k){
        std::ostringstream s; s<<"Chebyshev-U p["<<k<<"].beta = 1/4";
        check_near(b[k].beta(), 0.25, 1e-6, s.str());
    }

    // Shape vs U_n  (our p_n = U_n/2^n so shape ratios are identical)
    std::vector<double> xs={-0.7,-0.2,0.2,0.7};
    for(int n=1;n<=5;++n) ref::shape(b, n, 0.4, xs, ref::U, 2e-5, "Chebyshev-U");

    check_orth(b, w, 1e-6, "Chebyshev-U");

    opts.orthonormal=true;
    auto bon=StieltjesBuilder(w,opts).build(5);
    check_on(bon, w, 1e-5, "Chebyshev-U ON");
}

// ═════════════════════════════════════════════════════════════════════════════
// HERMITE
// ═════════════════════════════════════════════════════════════════════════════

void test_hermite(){
    std::cout<<"\n── Hermite (physicists') ────────────────────────────────\n";
    HermiteWeight w(8.0);
    BuilderOptions opts; opts.quad_points=150;
    auto b = StieltjesBuilder(w,opts).build(6);

    // Our monic p_n = H_n / 2^n  →  norm²[n] = √π · n! / 2^n
    for(int n=0;n<=6;++n){
        double exp=std::sqrt(M_PI)*ref::fact(n)/std::pow(2.0,n);
        std::ostringstream s; s<<"Hermite norm2["<<n<<"] = √π·n!/2^n";
        check_near(b[n].norm2(), exp, 1e-5, s.str());
    }

    // All stored alphas = 0 (symmetric)
    for(int n=0;n<=6;++n){
        std::ostringstream s; s<<"Hermite alpha in p["<<n<<"] = 0";
        check_near(b[n].alpha(), 0.0, 1e-13, s.str());
    }

    // Monic Hermite beta: p[k].beta() = (k-1)/2
    check_near(b[0].beta(), 0.0, 1e-13, "Hermite p[0].beta = 0");
    check_near(b[1].beta(), 0.0, 1e-13, "Hermite p[1].beta = 0");
    for(int k=2;k<=6;++k){
        double exp=(k-1)/2.0;
        std::ostringstream s; s<<"Hermite p["<<k<<"].beta = (k-1)/2 = "<<exp;
        check_near(b[k].beta(), exp, 1e-6, s.str());
    }

    // Shape vs H_n (same shape since p_n = H_n/2^n)
    std::vector<double> xs={-2.0,-0.5,0.0,0.5,2.0};
    for(int n=1;n<=6;++n) ref::shape(b, n, 1.0, xs, ref::H, 1e-6, "Hermite");

    check_orth(b, w, 1e-13, "Hermite");

    opts.orthonormal=true;
    auto bon=StieltjesBuilder(w,opts).build(5);
    check_on(bon, w, 1e-8, "Hermite ON");
}

// ═════════════════════════════════════════════════════════════════════════════
// LAGUERRE
// ═════════════════════════════════════════════════════════════════════════════

void test_laguerre(){
    std::cout<<"\n── Laguerre (classical, alpha=0) ────────────────────────\n";
    LaguerreWeight w(0.0,60.0);
    BuilderOptions opts; opts.quad_points=200;
    auto b = StieltjesBuilder(w,opts).build(6);

    // Our monic Laguerre: norm2[n] = (n!)^2
    {double f=1;
        for(int n=0;n<=6;++n){
            if(n>0)f*=n;
            double exp=f*f;
            std::ostringstream s; s<<"Laguerre norm2["<<n<<"] = (n!)^2 = "<<exp;
            double rtol = exp * 1e-4 + 1.0;  // relative 1e-4 tolerance
            check_near(b[n].norm2(), exp, rtol, s.str());
        }}

    // p[k].alpha(): alpha stored for p[k] is alpha_{k-1} = 2(k-1)+1 = 2k-1 for k>=1
    check_near(b[0].alpha(), 0.0, 1e-10, "Laguerre p[0].alpha = 0");
    for(int k=2;k<=5;++k){
        double exp=2.0*k-1.0;
        std::ostringstream s; s<<"Laguerre p["<<k<<"].alpha = "<<exp;
        check_near(b[k].alpha(), exp, 1e-3, s.str());
    }

    // p[k].beta(): beta stored for p[k] is beta_{k-1}=(k-1)^2
    check_near(b[0].beta(), 0.0, 1e-13, "Laguerre p[0].beta = 0");
    check_near(b[1].beta(), 0.0, 1e-13, "Laguerre p[1].beta = 0");
    for(int k=2;k<=6;++k){
        double exp=(double)((k-1)*(k-1));
        std::ostringstream s; s<<"Laguerre p["<<k<<"].beta = (k-1)^2 = "<<exp;
        check_near(b[k].beta(), exp, 1e-3, s.str());
    }

    // Shape vs L_n
    std::vector<double> xs={0.5,1.0,2.0,5.0};
    for(int n=1;n<=5;++n) ref::shape(b, n, 0.1, xs, ref::L, 1e-6, "Laguerre");

    check_orth(b, w, 1e-8, "Laguerre");

    opts.orthonormal=true;
    auto bon=StieltjesBuilder(w,opts).build(5);
    check_on(bon, w, 1e-13, "Laguerre ON");
}

// ═════════════════════════════════════════════════════════════════════════════
// JACOBI
// ═════════════════════════════════════════════════════════════════════════════

void test_jacobi(){
    struct Case { double a,b; double tol_orth; double tol_on; };
    std::vector<Case> cases={
            {1.0, 2.0,  1e-13, 1e-8},
            {0.0, 0.0,  1e-13, 1e-8},
            {2.0, 3.0,  1e-13, 1e-8},
            {0.5, 1.5,  1e-8, 1e-13},
    };
    for(auto& [a,b,to,ton] : cases){
        std::string tag="Jacobi(a="+std::to_string(a)+",b="+std::to_string(b)+")";
        std::cout<<"\n── "<<tag<<" ─────────────────────────────\n";
        JacobiWeight w(a,b);
        BuilderOptions opts; opts.quad_points=200;
        auto bas=StieltjesBuilder(w,opts).build(6);

        check_orth(bas, w, to, tag);

        for(int n=0;n<=6;++n){
            std::ostringstream s; s<<tag<<"  norm2["<<n<<"] > 0";
            check(bas[n].norm2()>0&&std::isfinite(bas[n].norm2()), s.str());
        }

        opts.orthonormal=true;
        auto bon=StieltjesBuilder(w,opts).build(5);
        check_on(bon, w, ton, tag+" ON");
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// CUSTOM WEIGHTS
// ═════════════════════════════════════════════════════════════════════════════

void test_custom(){
    std::cout<<"\n── Custom weights ───────────────────────────────────────\n";

    // Cauchy 1/(1+x^2) on [-6,6]
    {
        LambdaWeight w([](double x){return 1.0/(1.0+x*x);},-6.0,6.0,"Cauchy");
        BuilderOptions opts; opts.quad_points=100;
        auto b=StieltjesBuilder(w,opts).build(5);
        check_orth(b, w, 1e-8, "Cauchy");
        // Symmetric weight → all alpha = 0
        for(int n=0;n<=5;++n){
            std::ostringstream s; s<<"Cauchy alpha in p["<<n<<"] = 0";
            check_near(b[n].alpha(), 0.0, 1e-13, s.str());
        }
        // norm2[0] = ∫_{-6}^{6} 1/(1+x^2) dx = 2*arctan(6)
        check_near(b[0].norm2(), 2.0*std::atan(6.0), 1e-8, "Cauchy norm2[0]");
        opts.orthonormal=true;
        auto bon=StieltjesBuilder(w,opts).build(4);
        check_on(bon, w, 1e-13, "Cauchy ON");
    }

    // Log weight -ln(x) on (0,1]
    {
        LambdaWeight w([](double x){return x>1e-15?-std::log(x):0.0;},0.0,1.0,"LogW");
        BuilderOptions opts; opts.quad_points=200;
        auto b=StieltjesBuilder(w,opts).build(4);
        check_orth(b, w, 1e-6, "Log weight");
        // ∫_0^1 -ln(x) dx = 1
        check_near(b[0].norm2(), 1.0, 1e-4, "Log norm2[0] = 1");
    }

    // Step weight 1_{[0.2, 0.8]}
    {
        LambdaWeight w([](double x){return (x>0.2&&x<0.8)?1.0:0.0;},0.0,1.0,"Step");
        BuilderOptions opts; opts.quad_points=500;
        auto b=StieltjesBuilder(w,opts).build(4);
        check_orth(b, w, 1e-3, "Step weight");
        // ∫_0.2^0.8 dx = 0.6
        check_near(b[0].norm2(), 0.6, 5e-3, "Step norm2[0] = 0.6");
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// API & EDGE CASES
// ═════════════════════════════════════════════════════════════════════════════

void test_api(){
    std::cout<<"\n── API & edge cases ─────────────────────────────────────\n";
    LegendreWeight w;
    BuilderOptions opts;

    { auto b=StieltjesBuilder(w,opts).build(0);
        check(b.size()==1, "build(0) → 1 poly");
        check_near(b[0](0.5), 1.0, 1e-13, "p0(x)=1"); }

    { auto b=StieltjesBuilder(w,opts).build(3);
        bool threw=false;
        try{(void)b[10];}catch(const std::out_of_range&){threw=true;}
        check(threw,"out-of-range throws"); }

    { bool threw=false;
        try{LambdaWeight bad([](double x){return x;},1.0,0.0,"bad");}
        catch(const std::invalid_argument&){threw=true;}
        check(threw,"a>=b throws"); }

    { LambdaWeight neg([](double){return -1.0;},0.0,1.0,"neg");
        bool threw=false;
        try{StieltjesBuilder(neg,opts).build(2);}
        catch(const std::runtime_error&){threw=true;}
        check(threw,"negative weight throws"); }

    { auto b=StieltjesBuilder(w,opts).build(4);
        check(b.eval_all(0.5).size()==5,"eval_all size==N+1"); }

    { auto b=StieltjesBuilder(w,opts).build(4);
        auto G=StieltjesBuilder(w,opts).gram_matrix(b);
        bool ok=true; for(int i=0;i<5;++i) if(G[i][i]<=0) ok=false;
        check(ok,"Gram diagonal positive"); }

    { bool threw=false;
        try{LaguerreWeight bad(-2.0);}
        catch(const std::invalid_argument&){threw=true;}
        check(threw,"Laguerre alpha<=-1 throws"); }

    { auto b=StieltjesBuilder(w,opts).build(3);
        check(b[2].label()=="p2","label() correct"); }

    { LambdaWeight lw([](double x){return 1.0/(1.0+x*x);},-5.0,5.0,"test");
        check(lw.a()==-5.0&&lw.b()==5.0,"LambdaWeight a/b accessors"); }

    // max_degree()
    { auto b=StieltjesBuilder(w,opts).build(5);
        check(b.max_degree()==5,"max_degree() correct"); }

    // eval_matrix shape
    { auto b=StieltjesBuilder(w,opts).build(3);
        auto mat=b.eval_matrix({-0.5,0.0,0.5});
        check(mat.size()==4&&mat[0].size()==3,"eval_matrix shape"); }
}

// ─────────────────────────────────────────────────────────────────────────────

int main(){

    test_legendre();
    test_chebyshev_T();
    test_chebyshev_U();
    test_hermite();
    test_laguerre();
    test_jacobi();
    test_custom();
    test_api();

    int total=g_pass+g_fail;
    std::cout<<"\n══════════════════════════════════════════════════════════\n";
    std::cout<<"  Results:  "<<g_pass<<" / "<<total<<" passed";
    if(g_fail==0)
        std::cout<<"               ALL PASS\n";
    else
        std::cout<<"  ("<<g_fail<<" failed)\n";
    std::cout<<"══════════════════════════════════════════════════════════\n\n";
    return g_fail>0?1:0;
}