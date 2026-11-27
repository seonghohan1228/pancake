#include "case_model.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <tuple>

namespace gui {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kSteadyPseudoStep = 2e-5;  // s, as in the deck steady runs

double wrap_degrees(double deg) {
    double out = std::fmod(deg, 360.0);
    if (out < 0.0) out += 360.0;
    return out;
}

bool near(double a, double b, double rel = 1e-9) {
    return std::abs(a - b) <= rel * std::max({1.0, std::abs(a), std::abs(b)});
}

}  // namespace

bool uses_coupled(const CaseModel& m, const FluidPair* pair) {
    if (m.dissolved_gas && m.gas_model == GasModel::Transport) return true;
    // Steady FBNS film energy (finite-volume form) needs a gas-free film (pure
    // oil, deck O-runs); with refrigerant in the oil a steady thermal case runs
    // COUPLED (deck X-runs). Transient thermal uses the simplified FBNS form.
    return m.thermal && m.analysis == AnalysisType::Steady && effective_content(m, pair) > 0.0;
}

std::string unsupported_reason(const CaseModel& m, const FluidPair* pair) {
    const bool coupled = uses_coupled(m, pair);
    if (coupled && pair && pair->loaded() && !pair->has_gas_phase) {
        if (!(m.dissolved_gas && m.gas_model == GasModel::Transport))
            return "A steady thermal case with refrigerant in the oil uses gas transport, which "
                   "needs gas-phase data; this pair has none.";
        return "This pair has no gas-phase data; gas transport is not available.";
    }
    if (m.analysis == AnalysisType::Transient && coupled)
        return "Gas transport is available for steady cases only; transient cases use the "
               "pressure-plateau gas model.";
    if (m.feed.turns_with_shaft && (coupled || m.analysis != AnalysisType::Transient))
        return "A feed that turns with the shaft is available only in transient plateau cases.";
    if (coupled && m.dissolved_gas && pair && pair->loaded() && !(pair->gas_viscosity > 0.0))
        return "Gas transport needs the free-gas viscosity of the refrigerant in the fluid pair.";
    if (coupled && m.thermal && m.dissolved_gas && pair && pair->loaded() &&
        !pair->has_caloric_data())
        return "Thermal analysis with gas release needs caloric data in the fluid pair "
               "(dissolved heat capacity, desorption energy).";
    return "";
}

namespace {

double revolution_period(const CaseModel& m) {
    return m.speed != 0.0 ? 2.0 * kPi / std::abs(m.speed) : 0.0;
}

// --- key groups: each lists the case inputs it depends on ------------------

auto geometry_inputs(const CaseModel& m) {
    return std::tie(m.bearing_diameter, m.bearing_length, m.radial_clearance);
}
auto tilt_inputs(const CaseModel& m) { return std::tie(m.tilt_angle, m.tilt_direction_deg); }
auto speed_inputs(const CaseModel& m) { return std::tie(m.speed); }
auto fluid_inputs(const CaseModel& m) {
    return std::tie(m.fluid_pair, m.saturated, m.content, m.ambient_pressure,
                    m.ambient_temperature, m.film_temperature, m.dissolved_gas, m.gas_model,
                    m.thermal, m.analysis);
}
auto feed_inputs(const CaseModel& m) {
    return std::tuple_cat(std::tie(m.feed, m.bearing_length, m.speed), fluid_inputs(m));
}
auto motion_inputs(const CaseModel& m) {
    return std::tie(m.analysis, m.position, m.load, m.load_direction_deg, m.load_table,
                    m.moving_mass, m.fixed_eccentricity_ratio, m.fixed_attitude_deg,
                    m.radial_clearance);
}
auto roughness_inputs(const CaseModel& m) { return std::tie(m.roughness_shaft, m.roughness_bush); }
auto thermal_inputs(const CaseModel& m) {
    return std::tuple_cat(std::tie(m.housing, m.cross_film_layers, m.n_z, m.bearing_diameter),
                          fluid_inputs(m));
}
auto mesh_inputs(const CaseModel& m) { return std::tie(m.mesh_level, m.n_theta, m.n_z); }
auto time_inputs(const CaseModel& m) {
    return std::tie(m.analysis, m.speed, m.steps_per_revolution, m.revolutions,
                    m.outputs_per_revolution);
}
auto convergence_inputs(const CaseModel& m) { return std::tie(m.convergence); }

// --- group writers ----------------------------------------------------------

void write_geometry(const CaseModel& m, SimulationConfig& cfg) {
    cfg.R = 0.5 * m.bearing_diameter;
    cfg.L = m.bearing_length;
    cfg.c = m.radial_clearance;
}

void write_tilt(const CaseModel& m, SimulationConfig& cfg) {
    const double direction = m.tilt_direction_deg * kPi / 180.0;
    cfg.tilt_slope_x = m.tilt_angle * std::cos(direction);
    cfg.tilt_slope_y = m.tilt_angle * std::sin(direction);
}

void write_fluid_state(const CaseModel& m, const FluidPair* pair, SimulationConfig& cfg) {
    // Isothermal film temperature: ambient, unless set (two-step thermal route).
    const double temperature = (!m.thermal && m.film_temperature > 0.0) ? m.film_temperature
                                                                         : m.ambient_temperature;
    cfg.bc_z_south_type = cfg.bc_z_north_type = BCType::DIRICHLET;
    cfg.bc_z_south_val = cfg.bc_z_north_val = m.ambient_pressure;
    cfg.bc_z_south_thermal = cfg.bc_z_north_thermal = ThermalInflowMode::CONSTANT;
    cfg.bc_z_south_temperature = cfg.bc_z_north_temperature = temperature;
    cfg.temperature_initial = temperature;
    cfg.temperature_reference = temperature;

    const bool transport = uses_coupled(m, pair);
    if (transport) {
        // Gas transport (gas on: equilibrium release) or refrigerant-oil film
        // energy without release (gas off), as the deck X-runs.
        cfg.cavitation_model = CavitationModel::COUPLED;
        cfg.fluid_property_model = m.dissolved_gas ? FluidPropertyModel::TWO_PHASE
                                                   : FluidPropertyModel::SINGLE_PHASE;
        cfg.coupled_phase_equilibrium = m.dissolved_gas;
        // FBNS starting guess for the steady COUPLED solve; the solver accepts it
        // only with the finite-volume film energy (see coupled_fbns_seed_allowed).
        cfg.coupled_fbns_initialization = m.analysis == AnalysisType::Steady && m.thermal;
        cfg.cavitation_threshold = CavitationThreshold::SCALAR_PCAV;
        cfg.mixture_viscosity_model = MixtureViscosityModel::DUKLER;
        cfg.cavitated_film_model = CavitatedFilmModel::FULL_FILM;  // STRIATED is JFO/FBNS only
    } else {
        // FBNS: dissolved gas holds the cavitated film at the saturation pressure
        // (LOCAL_PSAT plateau); without it the film cavitates at the oil vapour
        // pressure. The dissolved content stays set either way, as in the deck runs.
        cfg.cavitation_model = CavitationModel::FBNS;
        cfg.fluid_property_model = FluidPropertyModel::SINGLE_PHASE;
        cfg.cavitation_threshold = m.dissolved_gas ? CavitationThreshold::LOCAL_PSAT
                                                   : CavitationThreshold::SCALAR_PCAV;
        cfg.steady_gas_model = SteadyGasModel::FROZEN;
        // Fixed for every FBNS case (plan R3): no shear from the bubble; only
        // the liquid at fraction theta shears (solver key STRIATED).
        cfg.cavitated_film_model = CavitatedFilmModel::STRIATED;
    }
    cfg.gas_pressure_coupling = GasPressureCoupling::NONE;
    if (!pair || !pair->loaded()) return;

    // R410A has a real-fluid gas model the solver allows only with COUPLED; the
    // FBNS plateau does not use the species, so it keeps the propane label of the
    // deck runs.
    // FBNS does not use the gas species; R410A is allowed only with COUPLED, so
    // FBNS keeps the propane label of the deck runs.
    cfg.dissolved_gas_species = (!transport && (!pair->has_gas_phase ||
                                                pair->solver_species == DissolvedGasSpecies::R410A))
                                    ? DissolvedGasSpecies::PROPANE
                                    : pair->solver_species;
    const double content = effective_content(m, pair);
    cfg.solution_table_file = pair->table_file;
    // Loaded tables travel with the config, as after the solver's table load,
    // so validate() and saturation_pressure() work on the GUI side.
    cfg.solubility_curve_table = pair->solubility;
    cfg.kinematic_composition_table = pair->viscosity;
    cfg.solubility_model = SolubilityModel::CURVE_TABLE;
    cfg.solution_table_clamp = true;
    cfg.kinematic_viscosity_log_interpolation = true;
    cfg.dissolved_gas_initial = content;
    cfg.dissolved_gas_max = 1.0;
    cfg.dissolved_gas_henry_coeff = 0.0;
    cfg.bc_z_south_dissolved_gas = cfg.bc_z_north_dissolved_gas = content;
    cfg.p_cav = pair->oil_vapour_pressure;

    const double rho = pair->density(temperature);
    cfg.rho = rho;
    cfg.mu = pair->dynamic_viscosity(content, temperature).value;
    cfg.rho_cp = rho * pair->specific_heat;
    cfg.thermal_conductivity = pair->thermal_conductivity;
    if (transport && pair->gas_viscosity > 0.0) cfg.mu_gas = pair->gas_viscosity;
    if (transport && pair->has_caloric_data()) {
        cfg.dissolved_heat_capacity = pair->dissolved_heat_capacity;
        cfg.desorption_internal_energy = pair->desorption_energy;
        cfg.gas_thermal_conductivity = pair->gas_thermal_conductivity;
        if (pair->gas_heat_capacity_cv > 0.0) cfg.gas_heat_capacity_cv = pair->gas_heat_capacity_cv;
    }
    if (m.thermal || transport) {
        cfg.density_model = DensityModel::LINEAR_TEMPERATURE;
        cfg.density_temperature_intercept = pair->density_intercept;
        cfg.density_temperature_slope = pair->density_slope;
        cfg.liquid_viscosity_model = LiquidViscosityModel::KINEMATIC_TABLE;
        cfg.liquid_specific_heat = pair->specific_heat;
    } else {
        // Isothermal FBNS: properties fixed at the ambient state (deck convention).
        cfg.density_model = DensityModel::CONSTANT;
        cfg.liquid_viscosity_model = LiquidViscosityModel::CONSTANT;
    }
}

void write_feed(const CaseModel& m, const FluidPair* pair, SimulationConfig& cfg) {
    InletConfig inlet = cfg.inlets.empty() ? InletConfig{} : cfg.inlets.front();
    const FeedInput& feed = m.feed;
    inlet.theta = wrap_degrees(feed.angle_deg);
    if (feed.type == FeedType::Hole) {
        inlet.type = InletConfig::Type::CIRCULAR;
        inlet.size = 0.5 * feed.hole_diameter;
        inlet.z = feed.hole_axial_position >= 0.0 ? feed.hole_axial_position : 0.5 * m.bearing_length;
    } else {
        inlet.type = InletConfig::Type::GROOVE;
        inlet.size = feed.groove_width_deg;
        inlet.z = 0.0;
    }
    inlet.p_supply = m.ambient_pressure + feed.supply_overpressure;
    inlet.feeds_fresh_oil = true;
    inlet.t_supply = feed.supply_temperature > 0.0 ? feed.supply_temperature : m.ambient_temperature;
    if (pair && pair->loaded()) inlet.dissolved_gas = effective_content(m, pair);
    if (cfg.inlets.empty()) cfg.inlets.push_back(inlet);
    else cfg.inlets.front() = inlet;
    cfg.inlet_rotation_speed = feed.turns_with_shaft ? m.speed : 0.0;
}

void write_motion(const CaseModel& m, const SimulationConfig& base, bool new_case,
                  SimulationConfig& cfg) {
    // The solver applies the external load to the moving bearing; a load W on
    // the shaft is the same as -W on the bearing.
    const double direction = (m.load_direction_deg + 180.0) * kPi / 180.0;
    cfg.bearing_stiffness_x = cfg.bearing_stiffness_y = cfg.bearing_stiffness_z = 0.0;
    cfg.bearing_damping_x = cfg.bearing_damping_y = cfg.bearing_damping_z = 0.0;
    cfg.external_load_z = 0.0;
    cfg.bearing_initial_from_attitude = true;
    // Initial position: keep a valid one from the file (it is usually close),
    // else start at 0.8 c roughly opposite the load, as a Newton seed.
    const bool keep_position = !new_case && base.e > 0.0 && base.e < 0.95 * m.radial_clearance;
    if (!keep_position) {
        // Seed: journal displaced towards the shaft load, ~20 deg ahead.
        cfg.e = 0.8 * m.radial_clearance;
        cfg.attitude_angle_deg = wrap_degrees(m.load_direction_deg + 20.0);
    }
    if (m.analysis == AnalysisType::Steady) {
        cfg.solution_mode = SolutionMode::STEADY_STATE;
        cfg.load_model = LoadModel::CONSTANT;
        if (m.position == ShaftPosition::Fixed) {
            cfg.motion_model = MotionModel::STATIC;
            cfg.e = m.fixed_eccentricity_ratio * m.radial_clearance;
            cfg.attitude_angle_deg = m.fixed_attitude_deg;
            cfg.external_load_x = cfg.external_load_y = 0.0;
        } else {
            cfg.motion_model = MotionModel::STEADY_EQUILIBRIUM;
            cfg.external_load_x = m.load * std::cos(direction);
            cfg.external_load_y = m.load * std::sin(direction);
        }
        return;
    }
    cfg.solution_mode = SolutionMode::TRANSIENT;
    cfg.motion_model = MotionModel::MOVING_BEARING;
    cfg.motion_coupling = MotionCoupling::IMPLICIT;
    cfg.motion_time_method = TimeSteppingMethod::EULER_IMPLICIT;
    cfg.bearing_mass = m.moving_mass;
    if (m.load_table.file.empty()) {
        cfg.load_model = LoadModel::CONSTANT;
        cfg.external_load_x = m.load * std::cos(direction);
        cfg.external_load_y = m.load * std::sin(direction);
        return;
    }
    cfg.load_model = LoadModel::TABLE;
    cfg.load_table_file = m.load_table.file;
    cfg.load_table_abscissa = m.load_table.abscissa == LoadAbscissa::CrankAngle
                                  ? LoadTableAbscissa::PHASE_DEG
                                  : LoadTableAbscissa::TIME_S;
    cfg.load_table_format = m.load_table.format == LoadFormat::MagnitudeDirection
                                ? LoadTableFormat::POLAR
                                : LoadTableFormat::CARTESIAN;
    cfg.load_frequency = m.load_table.frequency_hz;
    cfg.load_phase_offset_deg = m.load_table.phase_offset_deg;
    cfg.load_table_scale = m.load_table.scale;
    cfg.external_load_x = cfg.external_load_y = 0.0;
}

int housing_axial_cells(int n_z) {
    for (int cells : {16, 12, 8, 6, 4, 3, 2, 1})
        if (n_z % cells == 0) return cells;
    return 1;
}

void write_thermal(const CaseModel& m, const FluidPair* pair, SimulationConfig& cfg) {
    const double temperature = m.ambient_temperature;
    cfg.thermal_ambient_temperature = temperature;
    cfg.shaft_ambient_temperature = temperature;
    cfg.bush_outer_temperature = temperature;
    if (!m.thermal) {
        cfg.temperature_model = TemperatureModel::ISOTHERMAL;
        cfg.cross_film_thermal_model = CrossFilmModel::NONE;
        cfg.bush_thermal_model = BushThermalModel::LOCAL;
        return;
    }
    if (m.analysis == AnalysisType::Transient) {
        // Simplified thermal model for transient cases (plan R2): the solver runs
        // the finite-volume film and 3-D housing in steady cases only.
        cfg.temperature_model = TemperatureModel::ENERGY_EQUATION;
        cfg.cross_film_thermal_model = CrossFilmModel::POLYNOMIAL;
        cfg.cross_film_profile_order = 4;
        cfg.cross_film_quadrature_points = 6;
        cfg.bush_thermal_model = BushThermalModel::CONDUCTION_2D;
        cfg.bush_outer_radius = m.housing.outer_diameter > 0.0 ? 0.5 * m.housing.outer_diameter : 0.0;
        cfg.bush_thermal_conductivity = m.housing.conductivity;
        cfg.bush_outer_boundary = BushOuterBoundary::CONVECTION;
        cfg.bush_outer_heat_transfer = m.housing.outer_heat_transfer;
        cfg.shaft_ambient_heat_transfer = m.housing.shaft_heat_transfer;
        cfg.shaft_ambient_area_ratio = m.housing.shaft_area_ratio;
        cfg.groove_mixing_model = GrooveMixingModel::PINNED;
        return;
    }
    // Film + housing + shaft, as in the deck thermal runs (O/X series).
    cfg.temperature_model = TemperatureModel::ENERGY_EQUATION;
    cfg.cross_film_thermal_model = CrossFilmModel::FINITE_VOLUME;
    cfg.cross_film_layers = m.cross_film_layers;
    cfg.cross_film_profile_order = 4;
    cfg.cross_film_quadrature_points = 6;
    cfg.cross_film_relaxation = 0.3;
    cfg.bush_thermal_model = BushThermalModel::CONDUCTION_3D;
    cfg.bush_outer_radius = m.housing.outer_diameter > 0.0 ? 0.5 * m.housing.outer_diameter : 0.0;
    cfg.bush_thermal_conductivity = m.housing.conductivity;
    cfg.bush_outer_boundary = BushOuterBoundary::CONVECTION;
    cfg.bush_outer_heat_transfer = m.housing.outer_heat_transfer;
    cfg.bush_radial_cells = 16;
    cfg.bush_axial_cells = housing_axial_cells(m.n_z);
    cfg.bush_end_natural_convection = true;
    cfg.shaft_ambient_heat_transfer = m.housing.shaft_heat_transfer;
    cfg.shaft_ambient_area_ratio = m.housing.shaft_area_ratio;
    // COUPLED film energy (deck X-runs): pinned feed temperature. Gas-free FBNS
    // film energy (deck O-runs): flux-balanced feed.
    cfg.groove_mixing_model = uses_coupled(m, pair) ? GrooveMixingModel::PINNED
                                                    : GrooveMixingModel::FLUX_BALANCED;
}

void write_mesh(const CaseModel& m, SimulationConfig& cfg) {
    int n_theta = m.n_theta, n_z = m.n_z;
    if (m.mesh_level != MeshLevel::Custom) mesh_level_cells(m.mesh_level, n_theta, n_z);
    cfg.n_theta_global = n_theta;
    cfg.n_z_global = n_z;
}

void write_time(const CaseModel& m, const SimulationConfig& base, bool new_case,
                SimulationConfig& cfg) {
    const double period = revolution_period(m);
    if (m.analysis == AnalysisType::Transient && period > 0.0) {
        cfg.dt = period / std::max(1, m.steps_per_revolution);
        cfg.end_t = m.revolutions * period;
        cfg.write_interval = period / std::max(1, m.outputs_per_revolution);
        return;
    }
    // Steady: dt is only the COUPLED pseudo-step; keep a valid one from the file.
    const bool keep = !new_case && base.solution_mode == SolutionMode::STEADY_STATE && base.dt > 0.0;
    cfg.dt = keep ? base.dt : kSteadyPseudoStep;
    cfg.end_t = cfg.dt;
    cfg.write_interval = cfg.dt;
}

void write_convergence(const CaseModel& m, SimulationConfig& cfg) {
    const ConvergenceInput& c = m.convergence;
    cfg.steady_force_tolerance = c.force_tolerance;
    cfg.coupling_tolerance = c.outer_tolerance;
    cfg.coupling_max_iters = c.outer_max_iterations;
    cfg.steady_output_interval = std::max(1, c.field_view_interval);
    cfg.fbns_tol = c.inner_tolerance;
    cfg.fbns_max_iters = c.inner_max_iterations;
    cfg.coupled_nonlinear_tolerance = c.transport_tolerance;
    cfg.motion_coupling_tolerance = c.motion_tolerance;
    cfg.motion_coupling_position_tolerance = c.motion_tolerance;
    cfg.motion_coupling_max_iters = c.motion_max_iterations;
}

}  // namespace

void mesh_level_cells(MeshLevel level, int& n_theta, int& n_z) {
    switch (level) {
        case MeshLevel::Coarse: n_theta = 180; n_z = 48; return;
        case MeshLevel::Medium: n_theta = 360; n_z = 96; return;
        case MeshLevel::Fine: n_theta = 720; n_z = 192; return;
        case MeshLevel::Custom: return;
    }
}

const char* mesh_level_name(MeshLevel level) {
    switch (level) {
        case MeshLevel::Coarse: return "coarse";
        case MeshLevel::Medium: return "medium";
        case MeshLevel::Fine: return "fine";
        case MeshLevel::Custom: return "custom";
    }
    return "custom";
}

double effective_content(const CaseModel& m, const FluidPair* pair) {
    if (!m.saturated || !pair || !pair->loaded()) return m.content;
    return pair->saturated_content(m.ambient_pressure, m.ambient_temperature).value;
}

SimulationConfig new_case_template() {
    SimulationConfig cfg;
    cfg.output_verbosity = OutputVerbosity::DETAILED;
    cfg.diagnostics_interval = 1;
    cfg.output_write_flat = true;
    cfg.output_write_3d = false;
    cfg.log_outer_iters = false;
    cfg.steady_output_interval = 2;
    cfg.output_dir = "results";
    cfg.filename_prefix = "solution";
    cfg.omega_ramp_time = 0.0;
    return cfg;
}

namespace {

/// The solver's rule for coupled_fbns_initialization (SimulationConfig::validate):
/// steady COUPLED, finite-volume film energy, constant oil vapour pressure and no
/// prescribed free gas at the ends or feeds.
bool coupled_fbns_seed_allowed(const SimulationConfig& cfg) {
    if (cfg.solution_mode != SolutionMode::STEADY_STATE || cfg.cavitation_model != CavitationModel::COUPLED ||
        cfg.cross_film_thermal_model != CrossFilmModel::FINITE_VOLUME || cfg.oil_vapor_temperature_coeff != 0)
        return false;
    if (cfg.bc_z_south_free_gas_mass_fraction != 0 || cfg.bc_z_north_free_gas_mass_fraction != 0) return false;
    for (const auto& inlet : cfg.inlets)
        if (inlet.free_gas_mass_fraction != 0) return false;
    return true;
}

}  // namespace

SimulationConfig to_config(const CaseModel& m, const FluidPair* pair,
                           const SimulationConfig& base,
                           const std::optional<CaseModel>& original) {
    SimulationConfig cfg = base;
    const bool new_case = !original.has_value();
    const auto changed = [&](auto inputs) {
        return new_case || inputs(*original) != inputs(m);
    };
    if (changed([](const CaseModel& x) { return geometry_inputs(x); })) write_geometry(m, cfg);
    if (changed([](const CaseModel& x) { return tilt_inputs(x); })) write_tilt(m, cfg);
    if (changed([](const CaseModel& x) { return speed_inputs(x); })) cfg.omega = m.speed;
    if (changed([](const CaseModel& x) { return fluid_inputs(x); })) write_fluid_state(m, pair, cfg);
    if (changed([](const CaseModel& x) { return feed_inputs(x); })) write_feed(m, pair, cfg);
    if (changed([](const CaseModel& x) { return motion_inputs(x); }))
        write_motion(m, base, new_case, cfg);
    if (changed([](const CaseModel& x) { return roughness_inputs(x); })) {
        cfg.roughness_rq_shaft = m.roughness_shaft;
        cfg.roughness_rq_bush = m.roughness_bush;
    }
    if (changed([](const CaseModel& x) { return thermal_inputs(x); })) write_thermal(m, pair, cfg);
    if (changed([](const CaseModel& x) { return mesh_inputs(x); })) write_mesh(m, cfg);
    if (changed([](const CaseModel& x) { return time_inputs(x); }))
        write_time(m, base, new_case, cfg);
    if (changed([](const CaseModel& x) { return convergence_inputs(x); }))
        write_convergence(m, cfg);
    // The seed is only a starting guess: drop it wherever the solver would reject
    // it (e.g. gas transport with thermal off), rather than failing the run.
    if (cfg.coupled_fbns_initialization && !coupled_fbns_seed_allowed(cfg)) cfg.coupled_fbns_initialization = false;
    return cfg;
}

GuiMeta read_gui_meta(const std::string& text) {
    GuiMeta meta;
    std::istringstream stream(text);
    std::string line;
    while (std::getline(stream, line)) {
        const auto start = line.find_first_not_of(" \t");
        if (start == std::string::npos || line.compare(start, 5, "#@gui") != 0) continue;
        const std::string body = line.substr(start + 5);
        const auto eq = body.find('=');
        if (eq == std::string::npos) continue;
        auto trim = [](std::string s) {
            const auto a = s.find_first_not_of(" \t\r");
            const auto b = s.find_last_not_of(" \t\r");
            return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
        };
        meta[trim(body.substr(0, eq))] = trim(body.substr(eq + 1));
    }
    return meta;
}

GuiMeta case_meta(const CaseModel& m) {
    GuiMeta meta;
    meta["fluid_pair"] = m.fluid_pair;
    meta["saturated"] = m.saturated ? "true" : "false";
    meta["mesh_level"] = mesh_level_name(m.mesh_level);
    meta["cores"] = std::to_string(std::max(1, m.cores));
    return meta;
}

CaseRead from_config(const SimulationConfig& cfg, const GuiMeta& meta, const FluidPair* pair) {
    CaseRead out;
    CaseModel& m = out.model;
    const auto meta_value = [&](const char* key) {
        const auto it = meta.find(key);
        return it == meta.end() ? std::string() : it->second;
    };

    m.bearing_diameter = 2.0 * cfg.R;
    m.bearing_length = cfg.L;
    m.radial_clearance = cfg.c;
    m.tilt_angle = std::hypot(cfg.tilt_slope_x, cfg.tilt_slope_y);
    m.tilt_direction_deg = m.tilt_angle > 0.0
                               ? wrap_degrees(std::atan2(cfg.tilt_slope_y, cfg.tilt_slope_x) * 180.0 / kPi)
                               : 0.0;
    m.speed = cfg.omega;
    m.ambient_pressure = cfg.bc_z_south_val;
    if (!near(cfg.bc_z_south_val, cfg.bc_z_north_val))
        out.notes.push_back("The two bearing ends have different pressures; the GUI shows the z = 0 end.");

    m.thermal = cfg.temperature_model == TemperatureModel::ENERGY_EQUATION;
    m.ambient_temperature = m.thermal ? cfg.thermal_ambient_temperature : cfg.temperature_initial;
    m.cross_film_layers = cfg.cross_film_layers;
    m.housing.outer_diameter = cfg.bush_outer_radius > 0.0 ? 2.0 * cfg.bush_outer_radius : 0.0;
    m.housing.conductivity = cfg.bush_thermal_conductivity;
    m.housing.outer_heat_transfer = cfg.bush_outer_heat_transfer;
    m.housing.shaft_heat_transfer = cfg.shaft_ambient_heat_transfer;
    m.housing.shaft_area_ratio = cfg.shaft_ambient_area_ratio;

    // Feed: the first inlet; further inlets stay in the file untouched.
    if (!cfg.inlets.empty()) {
        const InletConfig& inlet = cfg.inlets.front();
        m.feed.angle_deg = inlet.theta;
        if (inlet.type == InletConfig::Type::CIRCULAR) {
            m.feed.type = FeedType::Hole;
            m.feed.hole_diameter = 2.0 * inlet.size;
            m.feed.hole_axial_position = near(inlet.z, 0.5 * cfg.L) ? -1.0 : inlet.z;
        } else {
            m.feed.type = FeedType::Groove;
            m.feed.groove_width_deg = inlet.size;
        }
        m.feed.supply_overpressure = inlet.p_supply - m.ambient_pressure;
        m.feed.supply_temperature =
            inlet.feeds_fresh_oil && !near(inlet.t_supply, m.ambient_temperature) ? inlet.t_supply : -1.0;
        if (cfg.inlets.size() > 1)
            out.notes.push_back(std::to_string(cfg.inlets.size() - 1) +
                                " further feed(s) in this file are kept but not shown.");
    }
    m.feed.turns_with_shaft = cfg.inlet_rotation_speed != 0.0;

    m.analysis = cfg.solution_mode == SolutionMode::TRANSIENT ? AnalysisType::Transient
                                                              : AnalysisType::Steady;
    m.position = cfg.motion_model == MotionModel::STATIC ? ShaftPosition::Fixed
                                                         : ShaftPosition::LoadBalance;
    m.load = std::hypot(cfg.external_load_x, cfg.external_load_y);
    m.load_direction_deg =
        m.load > 0.0 ? wrap_degrees(std::atan2(-cfg.external_load_y, -cfg.external_load_x) * 180.0 / kPi)
                     : 270.0;
    m.fixed_eccentricity_ratio = cfg.c > 0.0 ? cfg.e / cfg.c : 0.0;
    m.fixed_attitude_deg = cfg.attitude_angle_deg;
    m.moving_mass = cfg.bearing_mass;
    if (cfg.load_model == LoadModel::TABLE) {
        m.load_table.file = cfg.load_table_file;
        m.load_table.abscissa = cfg.load_table_abscissa == LoadTableAbscissa::PHASE_DEG
                                    ? LoadAbscissa::CrankAngle
                                    : LoadAbscissa::Time;
        m.load_table.format = cfg.load_table_format == LoadTableFormat::POLAR
                                  ? LoadFormat::MagnitudeDirection
                                  : LoadFormat::Components;
        m.load_table.frequency_hz = cfg.load_frequency;
        m.load_table.phase_offset_deg = cfg.load_phase_offset_deg;
        m.load_table.scale = cfg.load_table_scale;
    }
    if (cfg.load_model == LoadModel::ROTATING)
        out.notes.push_back("This file uses a rotating load; the GUI keeps it but cannot edit it.");

    m.dissolved_gas = cfg.cavitation_model == CavitationModel::COUPLED
                          ? cfg.fluid_property_model == FluidPropertyModel::TWO_PHASE
                          : cfg.cavitation_threshold == CavitationThreshold::LOCAL_PSAT;
    m.gas_model = cfg.cavitation_model == CavitationModel::COUPLED ? GasModel::Transport
                                                                    : GasModel::Plateau;
    if (cfg.cavitation_model != CavitationModel::FBNS &&
        cfg.cavitation_model != CavitationModel::COUPLED)
        out.notes.push_back("This file uses a cavitation model the GUI does not offer; it is kept.");

    m.content = cfg.dissolved_gas_initial;
    m.fluid_pair = meta_value("fluid_pair");
    if (m.fluid_pair.empty() && pair) m.fluid_pair = pair->id;
    const std::string saturated = meta_value("saturated");
    if (!saturated.empty()) {
        m.saturated = saturated == "true";
    } else if (pair && pair->loaded()) {
        const double c_sat = pair->saturated_content(m.ambient_pressure, m.ambient_temperature).value;
        m.saturated = near(c_sat, m.content, 1e-6);
    } else {
        m.saturated = false;
    }
    m.roughness_shaft = cfg.roughness_rq_shaft;
    m.roughness_bush = cfg.roughness_rq_bush;

    m.n_theta = cfg.n_theta_global;
    m.n_z = cfg.n_z_global;
    m.mesh_level = MeshLevel::Custom;
    for (MeshLevel level : {MeshLevel::Coarse, MeshLevel::Medium, MeshLevel::Fine}) {
        int nt = 0, nz = 0;
        mesh_level_cells(level, nt, nz);
        if (nt == m.n_theta && nz == m.n_z) m.mesh_level = level;
    }

    const double period = revolution_period(m);
    if (m.analysis == AnalysisType::Transient && period > 0.0 && cfg.dt > 0.0) {
        m.steps_per_revolution = static_cast<int>(std::lround(period / cfg.dt));
        m.revolutions = cfg.end_t / period;
        m.outputs_per_revolution =
            cfg.write_interval > 0.0 ? std::max(1, static_cast<int>(std::lround(period / cfg.write_interval))) : 1;
    }

    ConvergenceInput& c = m.convergence;
    c.force_tolerance = cfg.steady_force_tolerance;
    c.outer_tolerance = cfg.coupling_tolerance;
    c.outer_max_iterations = cfg.coupling_max_iters;
    if (cfg.steady_output_interval > 0) c.field_view_interval = cfg.steady_output_interval;
    c.inner_tolerance = cfg.fbns_tol;
    c.inner_max_iterations = cfg.fbns_max_iters;
    c.transport_tolerance = cfg.coupled_nonlinear_tolerance;
    c.motion_tolerance = cfg.motion_coupling_tolerance;
    c.motion_max_iterations = cfg.motion_coupling_max_iters;

    const std::string cores = meta_value("cores");
    m.cores = cores.empty() ? 1 : std::max(1, std::atoi(cores.c_str()));
    const std::string level = meta_value("mesh_level");
    if (level == "custom") m.mesh_level = MeshLevel::Custom;

    // Settings the GUI does not show: compare with what a fresh GUI case would write.
    const SimulationConfig fresh = to_config(m, pair, new_case_template());
    const auto mine = config_key_values(cfg.to_config_text());
    const auto theirs = config_key_values(fresh.to_config_text());
    // File paths are GUI inputs (fluid pair, load table) that a saved case stores
    // relative to itself, so their text differs from a fresh case without meaning.
    static const char* const kIgnored[] = {"output_dir",          "filename_prefix", "e",
                                           "attitude_angle_deg",  "solution_table_file",
                                           "load_table_file"};
    for (const auto& [key, value] : mine) {
        if (std::find(std::begin(kIgnored), std::end(kIgnored), key) != std::end(kIgnored)) continue;
        const auto it = theirs.find(key);
        if (it == theirs.end() || it->second != value) out.hidden_keys.push_back(key);
    }
    return out;
}

std::map<std::string, std::string> config_key_values(const std::string& text) {
    std::map<std::string, std::string> values;
    std::istringstream stream(text);
    std::string line;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto start = line.find_first_not_of(" \t");
        if (start == std::string::npos || line[start] == '#') continue;
        const auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = line.substr(start, eq - start);
        key.erase(key.find_last_not_of(" \t") + 1);
        std::string value = line.substr(eq + 1);
        const auto a = value.find_first_not_of(" \t");
        value = a == std::string::npos ? std::string() : value.substr(a);
        value.erase(value.find_last_not_of(" \t") + 1);
        auto& slot = values[key];
        slot = slot.empty() ? value : slot + "\n" + value;
    }
    return values;
}

}  // namespace gui
