# Inject with -DCMAKE_PROJECT_INCLUDE=/absolute/path/to/this/file when building
# either historical or repaired pancake. This leaves both solver sources intact.
add_executable(phase_exchange_probe "${CMAKE_CURRENT_LIST_DIR}/phase_exchange_probe.cpp")
target_link_libraries(phase_exchange_probe PRIVATE pancake_core)
