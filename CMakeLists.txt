cmake_minimum_required(VERSION 3.20)
project(pancake_bubble)

# Configuration
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Dependencies
find_package(MPI REQUIRED)

# PETSc
set(PETSC_DIR "/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real"
    CACHE PATH "PETSc installation directory")
find_path(PETSC_INCLUDE_DIR petsc.h HINTS "${PETSC_DIR}/include" REQUIRED)
find_library(PETSC_LIBRARY petsc_real HINTS "${PETSC_DIR}/lib" REQUIRED)
message(STATUS "PETSc include: ${PETSC_INCLUDE_DIR}")
message(STATUS "PETSc library: ${PETSC_LIBRARY}")

set(BUBBLE_SOURCES
    src/bubble_config.cpp         # BubbleConfig: all solver parameters
    src/stereo_mesh.cpp           # Inverse stereographic mesh + geometry engine
    src/mask.cpp                  # Circular domain mask + ghost blanking
    src/communicator.cpp          # MPI ghost exchange: open (non-periodic) boundaries
    src/bubble_io.cpp             # VTK unstructured (vtu/pvtu/pvd) output
    src/field.cpp                 # Field/Fields: scalar field with ghost layers
    src/bubble_linear_system.cpp  # PETSc-backed 5-point sparse linear system
    src/surface_operators.cpp     # ddt, laplacian, divergence, add_source, gradient
    src/thickness_transport.cpp   # Gravity velocity, face fluxes, transport solve
    src/bubble_main.cpp           # Entry point and time loop
)
file(GLOB HEADERS "src/*.hpp")

add_executable(pancake_bubble ${BUBBLE_SOURCES} ${HEADERS})

target_include_directories(pancake_bubble PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${PETSC_INCLUDE_DIR})

target_link_libraries(pancake_bubble PRIVATE MPI::MPI_CXX ${PETSC_LIBRARY})

# Tests
add_subdirectory(tests)

# Info
message(STATUS "MPI C++ Compiler: ${MPI_CXX_COMPILER}")
