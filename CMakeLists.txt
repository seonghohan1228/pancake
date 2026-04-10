cmake_minimum_required(VERSION 3.20)
project(pancake)

# Configuration
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Dependencies
find_package(MPI REQUIRED)

# Sources
file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "src/*.hpp")

# Targets
add_executable(pancake ${SOURCES} ${HEADERS})

# Include src/ directory so local headers can be found
target_include_directories(pancake PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Link MPI libraries
target_link_libraries(pancake PRIVATE MPI::MPI_CXX)

# Info
message(STATUS "MPI C++ Compiler: ${MPI_CXX_COMPILER}")
