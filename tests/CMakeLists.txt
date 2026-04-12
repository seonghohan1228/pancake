set(PETSC_DIR "/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real"
    CACHE PATH "PETSc installation directory")
find_path(PETSC_INCLUDE_DIR petsc.h HINTS "${PETSC_DIR}/include" REQUIRED)
find_library(PETSC_LIBRARY petsc_real HINTS "${PETSC_DIR}/lib" REQUIRED)

# Shared sources (everything except main.cpp)
file(GLOB SRC_FILES "${CMAKE_SOURCE_DIR}/src/*.cpp")
list(FILTER SRC_FILES EXCLUDE REGEX ".*main\\.cpp$")

function(add_pancake_test name)
    add_executable(${name} ${name}.cpp ${SRC_FILES})
    target_include_directories(${name} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${PETSC_INCLUDE_DIR})
    target_link_libraries(${name} PRIVATE MPI::MPI_CXX ${PETSC_LIBRARY})
    target_compile_features(${name} PRIVATE cxx_std_20)
endfunction()

add_pancake_test(test_linear_system)
add_pancake_test(test_fvm)
add_pancake_test(test_eos)
add_pancake_test(test_elrod_a1)
add_pancake_test(test_elrod_a2)
add_pancake_test(test_inlets)
add_pancake_test(test_velocity)
add_pancake_test(test_bc)
