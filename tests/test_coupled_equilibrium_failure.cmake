file(MAKE_DIRECTORY "${TEST_DIR}")
set(config "R = 0.05
c = 0.000152
e = 0.00006
L = 0.08
n_theta_global = 24
n_z_global = 4
omega = 10
mu = 0.0277
rho = 860
p_cav = 100
bc_z_south_val = 1000000
bc_z_north_val = 1000000
inlet_groove = 90 30 1200000
solution_mode = STEADY_STATE
cavitation_model = COUPLED
mixture_viscosity_model = DUKLER
fluid_property_model = SINGLE_PHASE
motion_model = STEADY_EQUILIBRIUM
external_load_y = 200
coupling_max_iters = 1
temperature_model = ISOTHERMAL
output_dir = ${TEST_DIR}/results
")
file(WRITE "${TEST_DIR}/config.txt" "${config}")
execute_process(COMMAND "${SOLVER}" -c "${TEST_DIR}/config.txt"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE errors TIMEOUT 30)
if(NOT result EQUAL 4 OR NOT output MATCHES "Equilibrium: NOT CONVERGED")
    message(FATAL_ERROR "Capped coupled load must fail with code 4: ${result}\n${output}\n${errors}")
endif()
file(STRINGS "${TEST_DIR}/results/coupled_summary.csv" summary)
list(GET summary 1 row)
if(NOT row MATCHES "^0,1,0,")
    message(FATAL_ERROR "Coupled cap must report unaccepted, cap exhausted, no rollback: ${row}")
endif()
# Reuse the just-computed fluid force as an independent load at its exact position.
file(STRINGS "${TEST_DIR}/results/equilibrium_history.csv" history)
list(GET history 1 row)
string(REPLACE "," ";" columns "${row}")
list(GET columns 2 x)
list(GET columns 3 y)
list(GET columns 9 fluid_x)
list(GET columns 10 fluid_y)
string(REPLACE "coupling_max_iters = 1" "coupling_max_iters = 10" config "${config}")
foreach(axis x y)
    if(fluid_${axis} MATCHES "^-")
        string(SUBSTRING "${fluid_${axis}}" 1 -1 load_${axis})
    else()
        set(load_${axis} "-${fluid_${axis}}")
    endif()
endforeach()
string(REPLACE "external_load_y = 200" "external_load_y = ${load_y}\nexternal_load_x = ${load_x}" config "${config}")
string(APPEND config "bearing_initial_from_attitude = false\nbearing_initial_x = ${x}\nbearing_initial_y = ${y}\ncoupling_relaxation = 0.3\n")
file(WRITE "${TEST_DIR}/config.txt" "${config}")
execute_process(COMMAND "${SOLVER}" -c "${TEST_DIR}/config.txt"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE errors TIMEOUT 30)
if(NOT result EQUAL 0 OR NOT output MATCHES "Equilibrium: CONVERGED")
    message(FATAL_ERROR "Known coupled equilibrium must succeed: ${result}\n${output}\n${errors}")
endif()
file(STRINGS "${TEST_DIR}/results/coupled_summary.csv" summary)
list(GET summary 1 row)
if(NOT row MATCHES "^1,0,0,")
    message(FATAL_ERROR "Converged coupled run must report strict acceptance: ${row}")
endif()
