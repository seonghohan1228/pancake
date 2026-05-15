set(PETSC_DIR "/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real"
    CACHE PATH "PETSc installation directory")
find_path(PETSC_INCLUDE_DIR petsc.h HINTS "${PETSC_DIR}/include" REQUIRED)
find_library(PETSC_LIBRARY petsc_real HINTS "${PETSC_DIR}/lib" REQUIRED)

# Bubble solver shared sources (mirrors CMakeLists.txt; excludes entry point)
set(BUBBLE_SRC_FILES
    ${CMAKE_SOURCE_DIR}/src/bubble_config.cpp
    ${CMAKE_SOURCE_DIR}/src/stereo_mesh.cpp
    ${CMAKE_SOURCE_DIR}/src/mask.cpp
    ${CMAKE_SOURCE_DIR}/src/communicator.cpp
    ${CMAKE_SOURCE_DIR}/src/bubble_io.cpp
    ${CMAKE_SOURCE_DIR}/src/field.cpp
    ${CMAKE_SOURCE_DIR}/src/bubble_linear_system.cpp
    ${CMAKE_SOURCE_DIR}/src/surface_operators.cpp
    ${CMAKE_SOURCE_DIR}/src/thickness_transport.cpp
)

function(add_bubble_test name)
    add_executable(${name} ${name}.cpp ${BUBBLE_SRC_FILES})
    target_include_directories(${name} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${PETSC_INCLUDE_DIR})
    target_link_libraries(${name} PRIVATE MPI::MPI_CXX ${PETSC_LIBRARY})
    target_compile_features(${name} PRIVATE cxx_std_20)
endfunction()

# Bubble tests — added per phase:
# Phase 2: test_stereo_mesh, test_mask
add_bubble_test(test_stereo_mesh)
add_bubble_test(test_mask)

# Phase 3: test_surface_operators, test_thickness_transport
add_bubble_test(test_surface_operators)
add_bubble_test(test_thickness_transport)

# Phase 4: test_simple
# Phase 5: test_energy, test_evaporation, test_marangoni
# Phase 6: test_disjoining
