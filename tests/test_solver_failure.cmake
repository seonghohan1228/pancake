# End-to-end contract: a failed Newton solve must have a nonzero process result.
file(MAKE_DIRECTORY "${TEST_DIR}")
set(config "${TEST_DIR}/config.txt")
file(WRITE "${config}" "R = 0.02\nc = 0.0001\ne = 0.00007\nL = 0.04\nn_theta_global = 24\nn_z_global = 8\nomega = 200\np_cav = 100000\nbc_z_south_val = 100000\nbc_z_north_val = 100000\nsolution_mode = STEADY_STATE\ncavitation_model = FBNS\nfbns_max_iters = 1\nfbns_tol = 1e-14\ncoupling_max_iters = 2\ntemperature_model = ISOTHERMAL\noutput_dir = ${TEST_DIR}/results\n")
execute_process(COMMAND "${SOLVER}" -c "${config}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE errors TIMEOUT 30)
if(NOT result EQUAL 4)
    message(FATAL_ERROR "Expected nonlinear failure exit 4, got ${result}: ${output}\n${errors}")
endif()
if(NOT output MATCHES "Simulation stopped" OR output MATCHES "Simulation finished successfully")
    message(FATAL_ERROR "Failure was reported as success: ${output}\n${errors}")
endif()
