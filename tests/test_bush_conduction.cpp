// Steady 2-D (r, theta) conduction in the bush annulus (src/bush_conduction.cpp,
// PHYSICS.md 13.1.2).
//
// The module is pure serial arithmetic on the global theta grid, so every rank
// runs the identical checks; that is also the guarantee the MPI caller relies on
// (all-reduce the coupling, then let every rank solve the same problem).
//
// What each test protects:
//   1. a purely radial problem is reproduced EXACTLY - the logarithmic face
//      conductance is the reason a 16-cell radial mesh is enough;
//   2. a theta-uniform problem reproduces the LOCAL 1-D closure exactly, so the
//      new model is a generalisation of the old one and not a different model;
//   3. a symmetric bore flux gives a symmetric bore temperature (no directional
//      bias leaked in from the sweep ordering or the periodic wrap);
//   4. a localised hot patch is spread circumferentially: lower peak, wider
//      support than the LOCAL closure on the SAME input. This is the physics the
//      wall-BC study found missing;
//   5. with an adiabatic outer surface the bore integral of the flux is zero
//      while the local flux is not - conduction only moves the heat around.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <mpi.h>

#include "bush_conduction.hpp"
#include "config.hpp"

static void check(bool cond, const char* msg) {
    if (!cond) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << msg << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

static constexpr double kRInner = 0.05;

static SimulationConfig bush_config() {
    SimulationConfig cfg;
    cfg.R = kRInner;
    cfg.bush_thermal_conductivity = 250.0;   // Ferron's bronze bush
    cfg.bush_outer_radius = 0.10;
    cfg.bush_outer_heat_transfer = 80.0;     // Ferron's outer free convection
    cfg.thermal_ambient_temperature = 313.15;
    cfg.bush_radial_cells = 16;
    cfg.bush_outer_boundary = BushOuterBoundary::CONVECTION;
    cfg.bush_thermal_model = BushThermalModel::CONDUCTION_2D;
    return cfg;
}

// A theta-independent bore flux makes the annulus one-dimensional, and the exact
// solution is the classical logarithmic shell in series with the outer film:
//   T(R) = T_amb + q R ln(R2/R)/k + q R/(R2 h).
// The discretisation uses the exact shell conductance on every radial face, so
// this must hold to round-off at ANY radial resolution -- including 1 cell.
static void test_uniform_flux_matches_the_analytic_annulus() {
    for (int n_r : {1, 4, 16, 64}) {
        for (bool fixed_outer : {false, true}) {
            SimulationConfig cfg = bush_config();
            cfg.bush_radial_cells = n_r;
            cfg.bush_outer_boundary = fixed_outer ? BushOuterBoundary::FIXED_TEMPERATURE
                                                  : BushOuterBoundary::CONVECTION;
            const int n_theta = 120;
            const double q = 4.0e4;  // W/m^2 into the bush, uniform

            BushConduction::Coupling coupling;
            coupling.a0.assign(n_theta, q);
            coupling.a1.assign(n_theta, 0.0);
            coupling.slope.assign(n_theta, 0.0);

            const BushConduction::Solution s =
                BushConduction::solve(coupling, n_theta, kRInner, cfg);
            check(s.converged, "the uniform-flux solve must converge");

            double expected = cfg.thermal_ambient_temperature +
                              q * kRInner * std::log(cfg.bush_outer_radius / kRInner) /
                                  cfg.bush_thermal_conductivity;
            if (!fixed_outer) {
                expected += q * kRInner /
                            (cfg.bush_outer_radius * cfg.bush_outer_heat_transfer);
            }
            double err = 0.0;
            for (int i = 0; i < n_theta; ++i) err = std::max(err, std::abs(s.alpha[i] - expected));
            // Round-off, not discretisation: measured 3e-13 relative at every
            // radial resolution, which is the point of the logarithmic face
            // conductance. An absolute 1e-10 on a ~570 K answer would be a
            // machine-precision assertion, so the bound is relative.
            check(err < 1.0e-11 * expected,
                  "a uniform bore flux must reproduce the analytic annulus temperature");
            check(std::abs(expected - cfg.thermal_ambient_temperature) > 1.0,
                  "the analytic check must be run on a non-trivial temperature rise");
        }
    }
}

// The LOCAL closure this replaces is q = G (T_b - T_amb) with G the series
// resistance of the wall and the outer film. When nothing varies with theta the
// 2-D solve can do nothing else, so it must land on the LOCAL answer exactly.
// This is the test that the Robin bore coupling and the half-cell elimination
// are wired up correctly, not just the interior operator.
static void test_theta_uniform_reduces_to_the_local_closure() {
    SimulationConfig cfg = bush_config();
    const int n_theta = 96;
    const double a0 = 3.1e6, slope = 6.0e3;  // q = a0 - slope * T_b

    BushConduction::Coupling coupling;
    coupling.a0.assign(n_theta, a0);
    coupling.a1.assign(n_theta, 0.0);
    coupling.slope.assign(n_theta, slope);
    const BushConduction::Solution s = BushConduction::solve(coupling, n_theta, kRInner, cfg);

    const double g = 1.0 / (kRInner * std::log(cfg.bush_outer_radius / kRInner) /
                                cfg.bush_thermal_conductivity +
                            kRInner / (cfg.bush_outer_radius * cfg.bush_outer_heat_transfer));
    const double expected = (a0 + g * cfg.thermal_ambient_temperature) / (slope + g);
    double err = 0.0;
    for (int i = 0; i < n_theta; ++i) err = std::max(err, std::abs(s.alpha[i] - expected));
    check(err < 1.0e-9 * expected,
          "a theta-uniform coupling must reproduce the LOCAL 1-D bush closure exactly");
    check(expected > cfg.thermal_ambient_temperature + 5.0,
          "the reduction check must be run on a non-trivial temperature rise");
}

// A bore flux symmetric under i -> n-1-i (a mirror about theta = 0) must give a
// symmetric bore temperature. Line relaxation sweeps theta in one direction and
// the periodic wrap is handled by a rank-one correction; either could smear a
// direction into the answer, and this is what would catch it.
static void test_symmetric_flux_gives_symmetric_temperature() {
    SimulationConfig cfg = bush_config();
    const int n_theta = 128;
    BushConduction::Coupling coupling;
    coupling.a0.assign(n_theta, 0.0);
    coupling.a1.assign(n_theta, 0.0);
    coupling.slope.assign(n_theta, 4.0e3);
    for (int i = 0; i < n_theta; ++i) {
        const double x = 2.0 * M_PI * (i + 0.5) / n_theta;
        coupling.a0[i] = 2.0e6 + 8.0e5 * std::cos(x) + 3.0e5 * std::cos(3.0 * x);
    }
    const BushConduction::Solution s = BushConduction::solve(coupling, n_theta, kRInner, cfg);

    double err = 0.0, spread = 0.0;
    double lo = s.alpha[0], hi = s.alpha[0];
    for (int i = 0; i < n_theta; ++i) {
        err = std::max(err, std::abs(s.alpha[i] - s.alpha[n_theta - 1 - i]));
        lo = std::min(lo, s.alpha[i]);
        hi = std::max(hi, s.alpha[i]);
    }
    spread = hi - lo;
    check(err < 1.0e-10 * (1.0 + spread), "a symmetric bore flux must give a symmetric bore temperature");
    check(spread > 1.0, "the symmetry check must be run on a genuinely varying field");
}

// The acid test in miniature. A localised hot patch on the bore is spread by
// conduction in the bush: against the LOCAL closure fed with the SAME coupling,
// the 2-D peak must be lower and the arc over which the temperature is elevated
// must be wider. A pure conduction operator is symmetric, so a symmetric patch
// broadens without moving; any circumferential SHIFT of the peak in a bearing
// comes from the asymmetry of the film flux, not from the bush.
static void test_hot_patch_is_spread_circumferentially() {
    SimulationConfig cfg = bush_config();
    const int n_theta = 360;
    const double base = 2.0e5, hot = 1.2e6, slope = 5.0e3;
    BushConduction::Coupling coupling;
    coupling.a0.assign(n_theta, base);
    coupling.a1.assign(n_theta, 0.0);
    coupling.slope.assign(n_theta, slope);
    for (int i = 150; i < 170; ++i) coupling.a0[i] = hot;  // 20 deg patch

    const BushConduction::Solution s = BushConduction::solve(coupling, n_theta, kRInner, cfg);

    const double g = 1.0 / (kRInner * std::log(cfg.bush_outer_radius / kRInner) /
                                cfg.bush_thermal_conductivity +
                            kRInner / (cfg.bush_outer_radius * cfg.bush_outer_heat_transfer));
    std::vector<double> local(n_theta, 0.0);
    for (int i = 0; i < n_theta; ++i) {
        local[i] = (coupling.a0[i] + g * cfg.thermal_ambient_temperature) / (slope + g);
    }

    // Effective support width (participation ratio) of the elevation above the
    // coldest point: (sum d)^2 / sum d^2, in cells. It returns exactly the patch
    // width for the LOCAL step and grows as the profile spreads, which the
    // half-maximum width does not once the peak itself has been cut.
    auto peak_and_width = [n_theta](const std::vector<double>& t, int& arg) {
        double lo = t[0], hi = t[0];
        arg = 0;
        for (int i = 0; i < n_theta; ++i) {
            lo = std::min(lo, t[i]);
            if (t[i] > hi) { hi = t[i]; arg = i; }
        }
        double sum = 0.0, sum_sq = 0.0;
        for (int i = 0; i < n_theta; ++i) {
            const double d = t[i] - lo;
            sum += d;
            sum_sq += d * d;
        }
        const double width = sum_sq > 0.0 ? sum * sum / sum_sq : 0.0;
        return std::pair<double, double>(hi - lo, width);
    };

    int arg_2d = 0, arg_local = 0;
    const auto [rise_2d, width_2d] = peak_and_width(s.alpha, arg_2d);
    const auto [rise_local, width_local] = peak_and_width(local, arg_local);

    check(rise_local > 1.0 && rise_2d > 0.0, "the patch must produce a measurable temperature rise");
    check(rise_2d < 0.9 * rise_local,
          "2-D bush conduction must cut the peak of a localised hot patch");
    check(width_2d > 2.0 * width_local,
          "2-D bush conduction must widen the arc over which the bush is hot");
    check(arg_2d >= 150 && arg_2d < 170,
          "a symmetric patch must keep its peak over the patch; conduction broadens, "
          "it does not advect");
    // The local closure sees the patch edge as a step; the 2-D bush cannot.
    double max_jump_2d = 0.0, max_jump_local = 0.0;
    for (int i = 0; i < n_theta; ++i) {
        const int ip = (i + 1) % n_theta;
        max_jump_2d = std::max(max_jump_2d, std::abs(s.alpha[ip] - s.alpha[i]));
        max_jump_local = std::max(max_jump_local, std::abs(local[ip] - local[i]));
    }
    check(max_jump_2d < 0.2 * max_jump_local,
          "the 2-D bore temperature must be smooth where the flux is discontinuous");
}

// An adiabatic outer surface leaves the bush a pure redistributor: the bore
// integral of the flux is zero to round-off, while the local flux is not. This is
// the discrete statement of global conservation for the assembled operator.
static void test_adiabatic_outer_conserves_energy() {
    SimulationConfig cfg = bush_config();
    cfg.bush_outer_heat_transfer = 0.0;  // adiabatic outer surface
    const int n_theta = 180;
    BushConduction::Coupling coupling;
    coupling.a1.assign(n_theta, 0.0);
    coupling.a0.assign(n_theta, 0.0);
    coupling.slope.assign(n_theta, 5.0e3);
    for (int i = 0; i < n_theta; ++i) {
        const double x = 2.0 * M_PI * (i + 0.5) / n_theta;
        coupling.a0[i] = 2.0e6 * (1.0 + 0.8 * std::sin(x) + 0.3 * std::sin(2.0 * x));
    }
    const BushConduction::Solution s = BushConduction::solve(coupling, n_theta, kRInner, cfg);

    double net = 0.0, worst = 0.0;
    for (int i = 0; i < n_theta; ++i) {
        const double q = coupling.a0[i] - coupling.slope[i] * s.alpha[i];
        net += q;
        worst = std::max(worst, std::abs(q));
    }
    check(worst > 1.0e4, "the conservation check needs a non-trivial local flux");
    check(std::abs(net) < 1.0e-8 * worst * n_theta,
          "an adiabatic bush must take no net heat while still exchanging locally");
}

// T_b = alpha + beta T_shaft must be the exact affine response of the linear
// problem: solving with the shaft contribution folded into a0 has to give the
// same answer as evaluating the two-part form.
static void test_shaft_slope_is_exact() {
    SimulationConfig cfg = bush_config();
    const int n_theta = 72;
    const double t_shaft = 331.7;
    BushConduction::Coupling split, folded;
    split.a0.assign(n_theta, 0.0);
    split.a1.assign(n_theta, 0.0);
    split.slope.assign(n_theta, 4.5e3);
    for (int i = 0; i < n_theta; ++i) {
        const double x = 2.0 * M_PI * (i + 0.5) / n_theta;
        split.a0[i] = 1.5e6 + 5.0e5 * std::cos(x);
        split.a1[i] = -2.0e3 - 4.0e2 * std::sin(x);
    }
    folded = split;
    for (int i = 0; i < n_theta; ++i) {
        folded.a0[i] = split.a0[i] + split.a1[i] * t_shaft;
        folded.a1[i] = 0.0;
    }

    const BushConduction::Solution a = BushConduction::solve(split, n_theta, kRInner, cfg);
    const BushConduction::Solution b = BushConduction::solve(folded, n_theta, kRInner, cfg);
    double err = 0.0, scale = 0.0;
    for (int i = 0; i < n_theta; ++i) {
        err = std::max(err, std::abs(a.alpha[i] + a.beta[i] * t_shaft - b.alpha[i]));
        scale = std::max(scale, std::abs(b.alpha[i]));
    }
    check(err < 1.0e-9 * scale,
          "alpha + beta T_shaft must reproduce the directly solved bore temperature");
    double max_beta = 0.0;
    for (int i = 0; i < n_theta; ++i) max_beta = std::max(max_beta, std::abs(a.beta[i]));
    check(max_beta > 1.0e-3, "the shaft slope must not be trivially zero in this setup");
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    test_uniform_flux_matches_the_analytic_annulus();
    test_theta_uniform_reduces_to_the_local_closure();
    test_symmetric_flux_gives_symmetric_temperature();
    test_hot_patch_is_spread_circumferentially();
    test_adiabatic_outer_conserves_energy();
    test_shaft_slope_is_exact();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "test_bush_conduction: all checks passed\n";

    MPI_Finalize();
    return 0;
}
