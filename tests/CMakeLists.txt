function(add_pancake_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE pancake_core)
    target_compile_features(${name} PRIVATE cxx_std_20)
    pancake_deploy_mingw_runtime(${name})
    if(name STREQUAL "test_eos" OR name STREQUAL "test_r410a")
        add_test(NAME ${name} COMMAND ${name})
    else()
        add_test(NAME ${name}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:${name}>)
    endif()
endfunction()

add_pancake_test(test_linear_system)
add_pancake_test(test_fvm)
add_pancake_test(test_eos)
add_pancake_test(test_elrod_a1)
add_pancake_test(test_elrod_a2)
add_test(NAME test_elrod_a2_mpi6
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6 $<TARGET_FILE:test_elrod_a2>)
add_test(NAME test_elrod_a2_mpi4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_elrod_a2>)
add_pancake_test(test_elrod_bc)
add_pancake_test(test_inlets)
add_pancake_test(test_velocity)
add_pancake_test(test_bc)
add_pancake_test(test_journal_motion)
add_pancake_test(test_forces)
add_pancake_test(test_contact)
add_pancake_test(test_energy)
# The steady energy equation is the operator the rank-dependent divergence hit
# first; exercise it on the odd split too.
add_test(NAME test_energy_mpi4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_energy>)
add_pancake_test(test_io)
add_pancake_test(test_fluid_properties)
add_pancake_test(test_gas_transport)
add_pancake_test(test_diagnostics)
add_pancake_test(test_mass_residual)
add_pancake_test(test_fbns)
add_test(NAME test_fbns_mpi6
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6 $<TARGET_FILE:test_fbns>)
# -n 1 and -n 4 are not redundant with -n 2 / -n 6. Preconditioner quality is a
# function of how the theta direction is cut, and the Elrod and Schur solves used
# to diverge at exactly these two counts while passing at 2 and 6. Leaving them
# unpoliced is how that regression survived several releases.
add_test(NAME test_fbns_mpi1
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_fbns>)
add_test(NAME test_fbns_mpi4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_fbns>)
add_pancake_test(test_schemes)
add_test(NAME test_schemes_mpi6
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6 $<TARGET_FILE:test_schemes>)
add_test(NAME test_schemes_mpi4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_schemes>)
add_pancake_test(test_property_tables)
target_compile_definitions(test_property_tables PRIVATE
    PANCAKE_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
add_pancake_test(test_bush_conduction)
add_pancake_test(test_cross_film)
add_pancake_test(test_resolved_cross_film)
foreach(ranks 1 4)
    add_test(NAME test_resolved_cross_film_mpi${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} $<TARGET_FILE:test_resolved_cross_film>)
endforeach()
# -n 4 covers the all-reduce/redistribute path of the 2-D bush on an odd split.
add_test(NAME test_cross_film_mpi4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_cross_film>)
add_pancake_test(test_config_keys)
target_compile_definitions(test_config_keys PRIVATE
    PANCAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
add_pancake_test(test_coupling)
add_pancake_test(test_gas_coupling)
add_pancake_test(test_physics_audit)
add_pancake_test(test_coupled_cavitation)
add_pancake_test(test_coupled_boundaries)
add_pancake_test(test_coupled_energy)
add_pancake_test(test_r410a)
add_executable(test_composition_viscosity test_composition_viscosity.cpp)
target_link_libraries(test_composition_viscosity PRIVATE pancake_core)
target_compile_definitions(test_composition_viscosity PRIVATE PANCAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
pancake_deploy_mingw_runtime(test_composition_viscosity)
add_test(NAME test_composition_viscosity COMMAND test_composition_viscosity)
add_executable(test_fw68d_properties test_fw68d_properties.cpp)
target_link_libraries(test_fw68d_properties PRIVATE pancake_core)
target_compile_definitions(test_fw68d_properties PRIVATE PANCAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
pancake_deploy_mingw_runtime(test_fw68d_properties)
add_test(NAME test_fw68d_properties COMMAND test_fw68d_properties)
foreach(ranks 1 4)
    add_test(NAME test_coupled_energy_mpi${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} $<TARGET_FILE:test_coupled_energy>)
    add_test(NAME test_coupled_cavitation_mpi${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} $<TARGET_FILE:test_coupled_cavitation>)
    add_test(NAME test_coupled_boundaries_mpi${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} $<TARGET_FILE:test_coupled_boundaries>)
endforeach()
add_test(NAME test_physics_audit_mpi1
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_physics_audit>)
add_test(NAME test_physics_audit_mpi4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_physics_audit>)
add_test(NAME test_solver_failure
    COMMAND ${CMAKE_COMMAND} -DSOLVER=$<TARGET_FILE:pancake>
        -DTEST_DIR=${CMAKE_CURRENT_BINARY_DIR}/solver_failure
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test_solver_failure.cmake)
add_test(NAME test_transient_time
    COMMAND ${CMAKE_COMMAND} -DSOLVER=$<TARGET_FILE:pancake>
        -DTEST_DIR=${CMAKE_CURRENT_BINARY_DIR}/transient_time
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test_transient_time.cmake)
# Aborts by design: passes only when the KSP guard prints its own diagnosis.
add_pancake_test(test_ksp_guard)
set_tests_properties(test_ksp_guard PROPERTIES
    PASS_REGULAR_EXPRESSION "linear solve 'ksp guard probe' failed")

add_pancake_test(test_steady_equilibrium)
foreach(ranks 1 4)
    add_test(NAME test_steady_equilibrium_mpi${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} $<TARGET_FILE:test_steady_equilibrium>)
endforeach()

add_test(NAME test_steady_equilibrium_failure
    COMMAND ${CMAKE_COMMAND} -DSOLVER=$<TARGET_FILE:pancake>
        -DTEST_DIR=${CMAKE_CURRENT_BINARY_DIR}/steady_equilibrium_failure
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test_steady_equilibrium_failure.cmake)
