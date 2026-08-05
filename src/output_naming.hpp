#pragma once

#include <array>
#include <string>

// Single source of truth for VTK output names. The internal field registry keeps
// descriptive names (pressure, temperature, velocity_theta, pressure_force_x, ...);
// the written VTK arrays use short OpenFOAM-style names and group Cartesian
// component triples into proper vector arrays. Both the writer (io.cpp) and the GUI
// preview reader use these helpers so they always agree.
namespace OutputNaming {

// Scalar internal-field name -> OpenFOAM-style output name. Anything not listed is
// written unchanged (rho, mu, h, heat_generation, ...). theta is the Elrod liquid
// fraction; it is exposed as film_content (could also be 'alpha' in OpenFOAM terms).
inline std::string scalar(const std::string& internal) {
    if (internal == "pressure")             return "p";
    if (internal == "temperature")          return "T";
    if (internal == "theta")                return "film_content";
    if (internal == "temperature_journal")  return "T_journal";
    if (internal == "temperature_bearing")  return "T_bearing";
    if (internal == "heat_flux_journal")    return "q_journal";
    if (internal == "heat_flux_bearing")    return "q_bearing";
    return internal;
}

// Optional film-surface fields of the cross-film thermal model: the two film
// surface temperatures T(eta=0), T(eta=1) [K] and the wall heat fluxes there
// [W/m^2], the latter positive INTO the fluid. They exist only when
// cross_film_thermal_model != NONE and the user selected them in output_fields.
// Listed here, next to their output names, so the config validator
// (SimulationConfig::validate), the module that writes them (CrossFilm) and the
// GUI reader all work from one set of selection keys.
inline const std::array<const char*, 4>& surface_fields() {
    static const std::array<const char*, 4> names = {{
        "temperature_journal",
        "temperature_bearing",
        "heat_flux_journal",
        "heat_flux_bearing",
    }};
    return names;
}

inline const char* velocity_name() { return "U"; }  // Cartesian velocity vector

// Cartesian vector groups assembled from per-component constant resultant fields.
// 'enable_*' are the config/output_field selection keys for the three components;
// the group is written if any of them is selected.
struct VectorGroup {
    const char* out;  // output vector array name
    const char* x;
    const char* y;
    const char* z;
};

inline const std::array<VectorGroup, 5>& vector_groups() {
    static const std::array<VectorGroup, 5> groups = {{
        {"Fp",   "pressure_force_x", "pressure_force_y", "pressure_force_z"},
        {"Fv",   "viscous_force_x",  "viscous_force_y",  "viscous_force_z"},
        {"F",    "fluid_force_x",    "fluid_force_y",    "fluid_force_z"},
        {"Fext", "external_load_x",  "external_load_y",  "external_load_z"},
        {"xB",   "bearing_x",        "bearing_y",        "bearing_z"},
    }};
    return groups;
}

// True if this internal field is a component consumed by a vector group (so the
// scalar writer skips it; velocity components are handled separately in io.cpp).
inline bool is_consumed_component(const std::string& internal) {
    for (const auto& g : vector_groups()) {
        if (internal == g.x || internal == g.y || internal == g.z) return true;
    }
    return false;
}

}  // namespace OutputNaming
