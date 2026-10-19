file(MAKE_DIRECTORY "${TEST_DIR}")
set(config "R = 0.0215
c = 0.0000425
e = 0.00003965941345966106
L = 0.04
attitude_angle_deg = 300
omega = 20.943951023931955
n_theta_global = 48
n_z_global = 8
mu = 0.008304798329665642
mu_gas = 0.00001265
rho = 935
p_cav = 0
bc_z_south_val = 1404224.696064497
bc_z_north_val = 1404224.696064497
inlet_groove = 90 13.298312524795712 1404224.696064497 300
inlet_species = 0 0.3 0
solution_mode = STEADY_STATE
cavitation_model = COUPLED
mixture_viscosity_model = DUKLER
fluid_property_model = SINGLE_PHASE
dissolved_gas_species = AIR
dissolved_gas_initial = 0.3
bc_z_south_dissolved_gas = 0.3
bc_z_north_dissolved_gas = 0.3
dissolved_gas_henry_coeff = 0.0000003
gas_mass_transfer_rate = 0
gas_absorption_rate = 0
temperature_model = ENERGY_EQUATION
temperature_initial = 300
temperature_reference = 300
property_reference_temperature = 300
liquid_specific_heat = 2010
dissolved_heat_capacity = 1600
desorption_internal_energy = 190000
gas_heat_capacity_cv = 718
gas_thermal_conductivity = 0.013
cross_film_thermal_model = FINITE_VOLUME
cross_film_layers = 4
bush_outer_heat_transfer = 100
shaft_ambient_heat_transfer = 1200
thermal_ambient_temperature = 300
groove_mixing_model = PINNED
coupling_max_iters = 2
coupling_relaxation = 0.3
coupling_tolerance = 1e-6
fbns_tol = 1e-10
coupled_fbns_initialization = true
log_outer_iters = true
")

function(reject name extra diagnosis)
    set(directory "${TEST_DIR}/${name}")
    file(MAKE_DIRECTORY "${directory}")
    file(WRITE "${directory}/config.txt" "${config}${extra}\noutput_dir = ${directory}/results\n")
    execute_process(COMMAND "${SOLVER}" -c "${directory}/config.txt"
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE errors TIMEOUT 15)
    if(result EQUAL 0 OR NOT "${output}${errors}" MATCHES "${diagnosis}")
        message(FATAL_ERROR "${name} must reject incompatible initialization: ${result}\n${output}\n${errors}")
    endif()
endfunction()
reject(transient "solution_mode = TRANSIENT" "coupled_fbns_initialization requires steady")
reject(non_coupled "cavitation_model = FBNS" "coupled_fbns_initialization requires steady")
reject(free_reservoir "bc_z_south_free_gas_mass_fraction = 0.01" "requires gas-free feeds and reservoirs")
reject(free_inlet "inlet_species = 0 0.3 0.01" "requires gas-free feeds and reservoirs")
reject(variable_vapor "oil_vapor_temperature_coeff = 0.01" "constant oil vapor pressure")
reject(invalid_boolean "coupled_fbns_initialization = maybe" "must be true or false")
reject(invalid_seed_tolerance "fbns_tol = 0" "fbns_tol must be positive")

foreach(branch off on)
    set(directory "${TEST_DIR}/${branch}")
    file(MAKE_DIRECTORY "${directory}")
    set(branch_config "${config}")
    if(branch STREQUAL "on")
        string(APPEND branch_config "fluid_property_model = TWO_PHASE\ncoupled_phase_equilibrium = true\n")
    endif()
    file(WRITE "${directory}/config.txt" "${branch_config}output_dir = ${directory}/results\n")
    execute_process(COMMAND "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" 4 "${SOLVER}" -c "${directory}/config.txt" -coupled_sub_pc_type lu
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE errors TIMEOUT 45)
    string(REGEX MATCHALL "Coupled FBNS initialization:" seeds "${output}")
    list(LENGTH seeds seed_count)
    if(NOT seed_count EQUAL 1 OR NOT (result EQUAL 0 OR result EQUAL 4))
        message(FATAL_ERROR "${branch} must seed once and reach final COUPLED checks: ${result}\n${output}\n${errors}")
    endif()
    file(STRINGS "${directory}/results/coupled_summary.csv" summary)
    list(GET summary 1 row)
    string(REPLACE "," ";" columns "${row}")
    list(GET columns 3 outer_iterations)
    list(GET columns 4 pressure_converged)
    list(GET columns 5 pressure_residual)
    if(NOT outer_iterations EQUAL 2 OR NOT pressure_converged EQUAL 1 OR NOT pressure_residual LESS 1e-10)
        message(FATAL_ERROR "${branch} must solve two unchanged COUPLED pressure sweeps: ${row}\n${output}\n${errors}")
    endif()
endforeach()
