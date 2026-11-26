#include "case_status.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace gui {

const char* node_name(Node node) {
    switch (node) {
        case Node::Geometry: return "Geometry";
        case Node::Fluid: return "Fluid";
        case Node::Operating: return "Operating conditions";
        case Node::Boundaries: return "Boundary conditions";
        case Node::Mesh: return "Mesh";
        case Node::Solution: return "Solution";
        case Node::Results: return "Results";
    }
    return "";
}

namespace {

void raise(NodeStatus& status, Light light, const std::string& reason) {
    if (static_cast<int>(light) > static_cast<int>(status.light)) {
        status.light = light;
        status.reasons.insert(status.reasons.begin(), reason);
    } else {
        status.reasons.push_back(reason);
    }
}

bool starts_with_any(const std::string& text, std::initializer_list<const char*> prefixes) {
    for (const char* prefix : prefixes)
        if (text.rfind(prefix, 0) == 0) return true;
    return false;
}

std::string mm(double metres) {
    char text[32];
    std::snprintf(text, sizeof(text), "%.3g mm", metres * 1e3);
    return text;
}

}  // namespace

Node node_for_message(const std::string& m) {
    if (starts_with_any(m, {"R ", "c ", "L ", "e ", "tilt", "inlet", "attitude"})) return Node::Geometry;
    if (starts_with_any(m, {"solution_table", "dissolved", "mu ", "mu_gas", "rho", "density",
                            "viscosity", "solubility", "p_cav", "Bubble", "Solution curve",
                            "kinematic", "liquid_viscosity", "liquid_specific", "desorption",
                            "gas_heat", "gas_thermal", "property", "The selected mixture"}))
        return Node::Fluid;
    if (starts_with_any(m, {"omega", "external_load", "load_", "bearing_mass", "STEADY_EQUILIBRIUM",
                            "motion_model", "bearing_initial"}))
        return Node::Operating;
    if (starts_with_any(m, {"bc_z", "bush_", "shaft_", "thermal_ambient", "Housing"}))
        return Node::Boundaries;
    if (starts_with_any(m, {"n_theta", "n_z", "cross_film_layers", "bush_axial"})) return Node::Mesh;
    return Node::Solution;
}

StatusArray compute_status(const StatusInputs& in) {
    StatusArray s;
    for (NodeStatus& node : s) node.light = Light::Green;
    const CaseModel& m = *in.model;
    auto& geometry = s[static_cast<int>(Node::Geometry)];
    auto& fluid = s[static_cast<int>(Node::Fluid)];
    auto& operating = s[static_cast<int>(Node::Operating)];
    auto& boundaries = s[static_cast<int>(Node::Boundaries)];
    auto& mesh = s[static_cast<int>(Node::Mesh)];
    auto& solution = s[static_cast<int>(Node::Solution)];

    // Geometry
    if (!(m.bearing_diameter > 0.0)) raise(geometry, Light::Yellow, "Enter the bearing diameter.");
    if (!(m.bearing_length > 0.0)) raise(geometry, Light::Yellow, "Enter the bearing length.");
    if (!(m.radial_clearance > 0.0)) raise(geometry, Light::Yellow, "Enter the radial clearance.");
    if (m.radial_clearance > 0.0 && m.bearing_diameter > 0.0 &&
        m.radial_clearance > 0.05 * m.bearing_diameter)
        raise(geometry, Light::Red, "The clearance is more than 5 % of the diameter; check the units.");
    if (m.feed.type == FeedType::Hole && m.feed.hole_diameter >= m.bearing_length)
        raise(geometry, Light::Red, "The feed hole is wider than the bearing length.");

    // Fluid
    if (m.fluid_pair.empty()) {
        raise(fluid, Light::Yellow, "Choose an oil-refrigerant pair.");
    } else if (!in.pair) {
        raise(fluid, Light::Red, "The fluid pair could not be loaded: " +
                                     (in.pair_error.empty() ? std::string("not found") : in.pair_error));
    } else if (in.estimates && !in.estimates->fluid_in_range) {
        raise(fluid, in.allow_out_of_range ? Light::Green : Light::Yellow,
              "The ambient state (content, temperature) is outside the fluid table; the solver "
              "would hold the table end values.");
    }
    if (!m.saturated && (m.content < 0.0 || m.content >= 1.0))
        raise(fluid, Light::Red, "The refrigerant content must be between 0 and 100 wt%.");

    // Operating conditions
    if (!(std::abs(m.speed) > 0.0)) raise(operating, Light::Yellow, "Enter the speed.");
    if (!(m.ambient_pressure > 0.0)) raise(operating, Light::Yellow, "Enter the ambient pressure.");
    if (!(m.ambient_temperature > 0.0)) raise(operating, Light::Yellow, "Enter the ambient temperature.");
    if (m.analysis == AnalysisType::Steady && m.position == ShaftPosition::LoadBalance &&
        !(m.load > 0.0))
        raise(operating, Light::Yellow, "Enter the load.");
    if (m.analysis == AnalysisType::Steady && m.position == ShaftPosition::Fixed &&
        !(m.fixed_eccentricity_ratio >= 0.0 && m.fixed_eccentricity_ratio < 1.0))
        raise(operating, Light::Red, "The eccentricity ratio must be between 0 and 1.");
    if (m.analysis == AnalysisType::Transient) {
        if (m.load_table.file.empty() && !(m.load > 0.0))
            raise(operating, Light::Yellow, "Import a load table or enter a constant load.");
        if (!m.load_table.file.empty() && !in.load_table_found)
            raise(operating, Light::Red, "The load table cannot be read: " + in.load_table_error);
        if (!m.load_table.file.empty() && m.load_table.abscissa == LoadAbscissa::CrankAngle &&
            !(m.load_table.frequency_hz > 0.0))
            raise(operating, Light::Yellow, "Enter the load cycle frequency.");
        if (!(m.moving_mass > 0.0)) raise(operating, Light::Yellow, "Enter the moving mass.");
        if (!(m.revolutions > 0.0)) raise(operating, Light::Yellow, "Enter the number of revolutions.");
    }

    // Boundary conditions
    if (m.thermal && m.housing.outer_diameter > 0.0 && m.housing.outer_diameter <= m.bearing_diameter)
        raise(boundaries, Light::Red,
              "The housing outer diameter must exceed the bearing diameter (" + mm(m.bearing_diameter) + ").");
    if (m.feed.supply_overpressure < 0.0)
        raise(boundaries, Light::Yellow, "The feed pressure is below ambient; check the feed zone.");

    // Mesh
    int n_theta = m.n_theta, n_z = m.n_z;
    if (m.mesh_level != MeshLevel::Custom) mesh_level_cells(m.mesh_level, n_theta, n_z);
    if (n_theta < 24 || n_z < 4) raise(mesh, Light::Red, "Use at least 24 x 4 cells.");
    if (m.feed.type == FeedType::Hole && m.bearing_diameter > 0.0 && m.bearing_length > 0.0) {
        const double cell_arc = 3.14159265358979323846 * m.bearing_diameter / std::max(1, n_theta);
        const double cell_axial = m.bearing_length / std::max(1, n_z);
        const double cells_across = m.feed.hole_diameter / std::max(cell_arc, cell_axial);
        if (cells_across < 4.0)
            raise(mesh, Light::Yellow, "The feed hole spans fewer than 4 cells; choose a finer mesh.");
    }

    // Solution: unsupported combinations, then the solver's own validation.
    const std::string unsupported = unsupported_reason(m, in.pair);
    if (!unsupported.empty()) raise(solution, Light::Red, unsupported);
    if (in.config && in.pair) {
        std::vector<std::string> errors, warnings;
        try {
            in.config->validate(errors, warnings);
        } catch (const std::exception& error) {
            errors.push_back(std::string("Solution curve table: ") + error.what());
        }
        if (unsupported.empty()) {
            for (const std::string& error : errors)
                raise(s[static_cast<int>(node_for_message(error))], Light::Red, error);
        }
    }

    s[static_cast<int>(Node::Results)].light = in.results_light;
    if (!in.results_reason.empty()) s[static_cast<int>(Node::Results)].reasons = {in.results_reason};
    return s;
}

bool ready_to_run(const StatusArray& status) {
    for (int i = 0; i < static_cast<int>(Node::Results); ++i)
        if (status[i].light != Light::Green) return false;
    return true;
}

std::string solver_form_text(const CaseModel& m, const FluidPair* pair) {
    std::string text;
    if (uses_coupled(m, pair)) {
        text = m.dissolved_gas ? "COUPLED gas transport (oil, dissolved and free refrigerant), "
                                 "equilibrium release"
                               : "COUPLED film solver without gas release";
    } else {
        text = m.dissolved_gas ? "FBNS (complementarity), cavitated film held at the saturation "
                                 "pressure (gas release plateau)"
                               : "FBNS (complementarity), vapour cavitation at the oil vapour pressure";
    }
    if (m.analysis == AnalysisType::Transient)
        text += "; transient, implicit load-position coupling (Newton)";
    else if (m.position == ShaftPosition::LoadBalance)
        text += "; steady load balance (Newton on the shaft position)";
    else
        text += "; steady, fixed shaft position";
    text += m.thermal ? "; film, housing and shaft temperatures (finite-volume film energy)"
                      : "; isothermal";
    return text;
}

}  // namespace gui
