cmake_minimum_required(VERSION 3.20)
project(pancake)

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

# Sources
file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "src/*.hpp")

# Targets
add_executable(pancake ${SOURCES} ${HEADERS})

# Include src/ directory and PETSc headers
target_include_directories(pancake PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${PETSC_INCLUDE_DIR})

# Link MPI and PETSc
target_link_libraries(pancake PRIVATE MPI::MPI_CXX ${PETSC_LIBRARY})

# Tests
add_subdirectory(tests)

# Info
message(STATUS "MPI C++ Compiler: ${MPI_CXX_COMPILER}")
