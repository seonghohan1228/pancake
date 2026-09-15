foreach(name stereo_mesh mask surface_operators thickness_transport)
    add_executable(test_bubble_${name} test_${name}.cpp)
    target_link_libraries(test_bubble_${name} PRIVATE pancake_bubble_core)
    pancake_deploy_mingw_runtime(test_bubble_${name})
    foreach(ranks 1 2 4)
        add_test(NAME test_bubble_${name}_mpi${ranks}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
                $<TARGET_FILE:test_bubble_${name}>)
        set_tests_properties(test_bubble_${name}_mpi${ranks} PROPERTIES TIMEOUT 60)
    endforeach()
endforeach()
add_test(NAME test_bubble_cli
    COMMAND ${CMAKE_COMMAND} -DSOLVER=$<TARGET_FILE:pancake_bubble>
        -DTEST_DIR=${CMAKE_CURRENT_BINARY_DIR}/cli
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test_cli.cmake)
