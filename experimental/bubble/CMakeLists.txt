add_library(pancake_bubble_core STATIC
    src/bubble_config.cpp
    src/stereo_mesh.cpp
    src/mask.cpp
    src/communicator.cpp
    src/bubble_io.cpp
    src/field.cpp
    src/bubble_linear_system.cpp
    src/surface_operators.cpp
    src/thickness_transport.cpp)
target_include_directories(pancake_bubble_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_features(pancake_bubble_core PUBLIC cxx_std_20)
target_link_libraries(pancake_bubble_core PUBLIC MPI::MPI_CXX)
if(MINGW)
    target_compile_definitions(pancake_bubble_core PUBLIC
        PETSC_FUNCTION_NAME_CXX=__func__
        PETSC_CXX_RESTRICT=__restrict__)
endif()
if(PETSC_TARGET)
    target_link_libraries(pancake_bubble_core PUBLIC ${PETSC_TARGET})
else()
    target_include_directories(pancake_bubble_core PUBLIC ${PETSC_INCLUDE_DIR})
    target_link_libraries(pancake_bubble_core PUBLIC ${PETSC_LIBRARY})
endif()
add_executable(pancake_bubble src/bubble_main.cpp)
target_link_libraries(pancake_bubble PRIVATE pancake_bubble_core)
pancake_deploy_mingw_runtime(pancake_bubble)
if(PANCAKE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
