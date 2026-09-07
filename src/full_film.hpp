#pragma once

#include <algorithm>

#include "config.hpp"
#include "field.hpp"

/// Single definition of the question "does this cell still hold an unbroken
/// liquid film?", which gates Poiseuille transport of momentum, enthalpy and
/// dissolved gas everywhere in the solver.
///
/// The test used to be spelled as the literal `theta >= 1.0` at fifteen sites in
/// `reynolds.cpp`, `energy.cpp`, `gas_transport.cpp` and `diagnostics.cpp`. That
/// literal is only correct while the liquid may fill the whole gap. With a
/// resolved free-gas phase the ceiling drops to `1 - alpha_gas`, and a ceiling
/// that has to move in fifteen places moves wrong: every switch in the domain
/// would turn off at once the instant any gas came out of solution, and the
/// solver would report a fully cavitated bearing with a converged residual.
/// Routing every user through `bound()` made that a one-line change, which is
/// the change `gas_volume_fraction` now carries.
/// See `docs/FBNS_GAS_COUPLING.md` sections 3 and 4.
///
/// Usage: prefer `is_full_film` / `is_full_film_face` for branches and
/// `gate` / `face_gate` where the caller multiplies a flux by a 1/0 switch.
/// `holds` is the value-level core, for callers that already have a film content
/// in hand (a saved iterate) or that need a tolerance band.
namespace FullFilm {

/// True when the selected models let a resolved free-gas phase displace liquid
/// from the gap, so the full-film ceiling drops below 1.
///
/// Only the FBNS complementarity path. The Elrod (JFO) path is deliberately left
/// out even under TWO_PHASE:
///
/// - With `gas_pressure_coupling = NONE` its free gas is volumetrically inert by
///   construction (PHYSICS.md 13.2.1). Case 05 had alpha_gas = 0 and does not
///   validate this gas-volume approximation.
/// - With `VOID_COUPLED` the Elrod operator already carries the moved ceiling in
///   `Reynolds`'s `VoidParams::theta_full`, where it feeds the barotropic EOS and
///   the switch function. Dropping the fifteen Poiseuille gates onto it as well
///   changes the formulation and requires separate verification. VOID_COUPLED
///   remains experimental and has no converged experimental validation.
///
/// Merging the two ceilings therefore stays open work. When it lands, this
/// predicate loses its `cavitation_model` term and `VoidParams::theta_full`
/// becomes a call to `bound`. See `docs/FBNS_GAS_COUPLING.md` section 5.
inline bool gas_displaces_liquid(const SimulationConfig& cfg) {
    return cfg.cavitation_model == CavitationModel::FBNS &&
           cfg.fluid_property_model == FluidPropertyModel::TWO_PHASE;
}

/// Volume fraction of the gap held by resolved free gas, as the full-film test
/// sees it. Zero whenever no gas phase displaces liquid, which makes every
/// non-FBNS path bit-identical to the pre-gas ceiling by construction.
///
/// Clamped to the same `[0, gas_alpha_max]` band `FluidProperties::update_gas_state`
/// writes and `Reynolds::void_params` reads, capped just below 1 so the ceiling
/// never reaches zero and strands the liquid balance without a film.
inline double gas_volume_fraction(const Fields& fields, const SimulationConfig& cfg,
                                  int i, int j) {
    if (!gas_displaces_liquid(cfg) || !fields.has("alpha_gas")) return 0.0;
    const double alpha_cap = std::min(std::max(cfg.gas_alpha_max, 0.0), 0.999);
    return std::clamp(fields["alpha_gas"](i, j), 0.0, alpha_cap);
}

/// Largest film content a cell can hold, i.e. the theta at and above which the
/// film is unbroken and the cell carries Poiseuille flow.
///
/// Exactly 1 unless a resolved free-gas phase occupies part of the gap, in which
/// case the liquid can only fill what is left: theta <= 1 - alpha_gas. Every
/// full-film test in the solver reads the ceiling from here, so moving it moves
/// all of them together -- which is the point, since moving it in only some
/// places reports a fully cavitated bearing with a converged residual
/// (docs/FBNS_GAS_COUPLING.md section 3).
///
/// alpha_gas is *lagged*: the operator-split gas path (`GasTransport::solve` then
/// `FluidProperties::update_gas_state`) runs after the flow solve, so within one
/// Newton solve this is a coefficient, not an unknown. That is what keeps
/// d(bound)/d(theta) = 0 and leaves the FBNS Jacobian blocks unchanged.
inline double bound(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return std::max(1.0 - gas_volume_fraction(fields, cfg, i, j), 0.0);
}

/// Value-level core test. `tol` widens the bound downwards for callers that must
/// not let linear-solver round-off flip the flag; it defaults to an exact test.
inline bool holds(double theta_value, double bound_value, double tol = 0.0) {
    return theta_value >= bound_value - tol;
}

/// True when cell (i,j) holds an unbroken liquid film.
inline bool is_full_film(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return holds(fields["theta"](i, j), bound(fields, cfg, i, j));
}

/// True when the face separating (i0,j0) from (i1,j1) is a full-film face: both
/// cells it joins must hold an unbroken film for the face to carry Poiseuille
/// flow, since either one ruptured cell breaks the pressure-driven path.
inline bool is_full_film_face(const Fields& fields, const SimulationConfig& cfg,
                              int i0, int j0, int i1, int j1) {
    return is_full_film(fields, cfg, i0, j0) && is_full_film(fields, cfg, i1, j1);
}

/// The cell test as the 1/0 switch `g` the flux expressions multiply through.
inline double gate(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    if (cfg.cavitation_model == CavitationModel::COUPLED) return 1.0;
    return is_full_film(fields, cfg, i, j) ? 1.0 : 0.0;
}

/// The face test as the 1/0 switch `g_face` the flux expressions multiply through.
inline double face_gate(const Fields& fields, const SimulationConfig& cfg,
                        int i0, int j0, int i1, int j1) {
    if (cfg.cavitation_model == CavitationModel::COUPLED) return 1.0;
    return is_full_film_face(fields, cfg, i0, j0, i1, j1) ? 1.0 : 0.0;
}

}  // namespace FullFilm
