#pragma once

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
/// Routing every user through `bound()` makes that a one-line change.
/// See `docs/FBNS_GAS_COUPLING.md` sections 3 and 4.
///
/// Usage: prefer `is_full_film` / `is_full_film_face` for branches and
/// `gate` / `face_gate` where the caller multiplies a flux by a 1/0 switch.
/// `holds` is the value-level core, for callers that already have a film content
/// in hand (a saved iterate) or that need a tolerance band.
namespace FullFilm {

/// Largest film content a cell can hold, i.e. the theta at and above which the
/// film is unbroken and the cell carries Poiseuille flow.
///
/// Exactly 1 for every model currently shipped: free gas is volumetrically inert
/// in the film-content equation, so the liquid may occupy the whole gap.
/// Admitting a resolved gas volume is the single line
///
///     return std::max(1.0 - gas_volume_fraction(fields, cfg, i, j), 0.0);
///
/// here, with no call-site edits. `Reynolds`'s `VoidParams::theta_full` is the
/// same ceiling as seen by the Elrod operator under `VOID_COUPLED`; the two
/// merge when that line lands.
inline double bound([[maybe_unused]] const Fields& fields,
                    [[maybe_unused]] const SimulationConfig& cfg,
                    [[maybe_unused]] int i, [[maybe_unused]] int j) {
    return 1.0;
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
    return is_full_film(fields, cfg, i, j) ? 1.0 : 0.0;
}

/// The face test as the 1/0 switch `g_face` the flux expressions multiply through.
inline double face_gate(const Fields& fields, const SimulationConfig& cfg,
                        int i0, int j0, int i1, int j1) {
    return is_full_film_face(fields, cfg, i0, j0, i1, j1) ? 1.0 : 0.0;
}

}  // namespace FullFilm
