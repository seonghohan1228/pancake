#include "r410a.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace R410A {
namespace {
// Lemmon (2003), DOI 10.1023/A:1025048800563; coefficient transcription
// checked against CoolProp 7.2.0 R410A.json. See docs/R410A_EOS.md.
struct Term { double n, t; int d, l; };
constexpr std::array<Term, 21> terms{{
    {.987252,.44,1,0}, {-1.03017,1.2,1,0}, {1.17666,2.97,1,0},
    {-.138991,2.95,2,0}, {.00302373,.2,5,0}, {-2.53639,1.93,1,1},
    {-1.9668,1.78,2,1}, {-.83048,3,3,1}, {.172477,.2,5,1},
    {-.261116,.74,5,1}, {-.0745473,3,5,1}, {.679757,2.1,1,2},
    {-.652431,4.3,1,2}, {.0553849,.25,4,2}, {-.071097,7,4,2},
    {-.000875332,4.7,9,2}, {.020076,13,2,3}, {-.0139761,16,2,3},
    {-.018511,25,4,3}, {.0171939,17,5,3}, {-.00482049,7.4,6,3}
}};
void check_temperature(double t) {
    if (!std::isfinite(t) || t < 200 || t > 500)
        throw std::domain_error("R410A EOS requires 200 <= T <= 500 K");
}
std::array<double, 2> evaluate(double rho, double t) {
    check_temperature(t);
    if (!std::isfinite(rho) || rho < 0)
        throw std::domain_error("R410A EOS requires nonnegative finite density");
    const double delta = rho / critical_density, tau = critical_temperature / t;
    thread_local double last_temperature = -1;
    thread_local std::array<double,21> temperature_terms{};
    if (last_temperature != t) {
        for (std::size_t i=0;i<terms.size();++i) temperature_terms[i]=terms[i].n*std::pow(tau,terms[i].t);
        last_temperature=t;
    }
    std::array<double,10> powers{}; powers[0]=1;
    for (std::size_t i=1;i<powers.size();++i) powers[i]=powers[i-1]*delta;
    const std::array<double,4> exponentials{{1,std::exp(-powers[1]),std::exp(-powers[2]),std::exp(-powers[3])}};
    double z = 1, derivative = 1;
    for (std::size_t i=0;i<terms.size();++i) {
        const auto& v=terms[i];
        const double dl = v.l ? powers[v.l] : 0;
        const double a = temperature_terms[i] * powers[v.d] * exponentials[v.l];
        const double b = v.d - v.l * dl;
        z += a * b;
        derivative += a * (b * (1 + b) - v.l * v.l * dl);
    }
    return {rho * gas_constant * t * z, gas_constant * t * derivative};
}
}
double dew_pressure(double t) {
    check_temperature(t);
    if (t >= critical_temperature) return std::numeric_limits<double>::infinity();
    const double x = 1 - t / critical_temperature;
    return 4901200 * std::exp(critical_temperature / t *
        (-7.4411*x + 1.9883*std::pow(x,1.6) - 2.4925*std::pow(x,2.4) - 3.2633*std::pow(x,5)));
}
double pressure(double rho, double t) { return evaluate(rho,t)[0]; }
double pressure_derivative(double rho, double t) { return evaluate(rho,t)[1]; }

namespace {
std::array<double,2> caloric(double rho, double t) {
    check_temperature(t);
    if (!std::isfinite(rho) || rho<0) throw std::domain_error("R410A caloric EOS requires finite nonnegative density");
    const double tau=critical_temperature/t, delta=rho/critical_density;
    // tau*alpha_tau and tau^2*alpha_tautau, including ideal and residual terms.
    double first=7.15807*tau-1+4.687575*std::pow(tau,-.1);
    double second=1-5.1563325*std::pow(tau,-.1);
    constexpr std::array<double,3> n{{2.0623,5.9751,1.5612}}, theta{{2.02326,5.00154,11.2484}};
    for (int k=0;k<3;++k) {
        const double x=theta[k]*tau, e=std::exp(-x), denominator=1-e;
        first+=n[k]*x*e/denominator;
        second-=n[k]*x*x*e/(denominator*denominator);
    }
    for (const auto& v : terms) {
        const double a=v.n*std::pow(tau,v.t)*std::pow(delta,v.d)*
            (v.l ? std::exp(-std::pow(delta,v.l)) : 1);
        first+=v.t*a; second+=v.t*(v.t-1)*a;
    }
    return {gas_constant*t*first,-gas_constant*second};
}
}
double internal_energy(double rho, double t) { return caloric(rho,t)[0]; }
double heat_capacity_cv(double rho, double t) { return caloric(rho,t)[1]; }
double gas_density(double p, double t) {
    check_temperature(t);
    if (!std::isfinite(p) || p < 0 || p > maximum_pressure)
        throw std::domain_error("R410A free-fluid EOS requires 0 <= p <= 50 MPa absolute");
    if (p > dew_pressure(t) * (1 + 1e-12))
        throw std::domain_error("R410A free phase is above its dew pressure; a refrigerant liquid/flash model is required");
    if (p == 0) return 0;
    // Follow the dilute branch and stop at the first mechanically unstable point.
    double lo = 0, hi = p / (gas_constant*t);
    for (int n = 0; ; ++n) {
        const auto v = evaluate(hi,t);
        if (v[0] >= p) break;
        if (n > 100 || v[1] <= 0)
            throw std::domain_error("R410A stable gas root could not be bracketed");
        hi *= 1.08;
    }
    double rho = (lo + hi)/2;
    for (int n = 0; n < 80; ++n) {
        const auto v = evaluate(rho,t);
        if (std::abs(v[0] - p) <= 2e-12 * p) return rho;
        if (v[0] > p) hi = rho; else lo = rho;
        const double next = rho - (v[0]-p)/v[1];
        rho = next > lo && next < hi ? next : (lo+hi)/2;
    }
    throw std::runtime_error("R410A gas density iteration failed");
}
}
