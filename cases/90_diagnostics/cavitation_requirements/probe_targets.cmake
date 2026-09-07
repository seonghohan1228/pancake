add_executable(kinetics_probe "${CMAKE_CURRENT_LIST_DIR}/kinetics_probe.cpp")
target_link_libraries(kinetics_probe PRIVATE pancake_core)
