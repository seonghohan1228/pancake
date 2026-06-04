cmake_minimum_required(VERSION 3.21)
project(pancake LANGUAGES CXX)

# Configuration
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(WIN32)
    set(PANCAKE_GUI_DEFAULT ON)
else()
    set(PANCAKE_GUI_DEFAULT OFF)
endif()

option(PANCAKE_BUILD_SOLVER "Build the MPI/PETSc solver and CLI executable" ON)
option(PANCAKE_BUILD_TESTS "Build solver tests" ON)
option(PANCAKE_BUILD_GUI "Build the native Windows config GUI" ${PANCAKE_GUI_DEFAULT})

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config.txt
    ${CMAKE_CURRENT_BINARY_DIR}/config.txt
    COPYONLY)

function(pancake_copy_default_config target_name)
    add_custom_command(TARGET ${target_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/config.txt"
            "$<TARGET_FILE_DIR:${target_name}>/config.txt"
        COMMENT "Copying config.txt next to target executable"
        VERBATIM)
endfunction()

function(pancake_deploy_mingw_runtime target_name)
    if(MINGW)
        get_filename_component(PANCAKE_MINGW_BIN "${CMAKE_CXX_COMPILER}" DIRECTORY)
        add_custom_command(TARGET ${target_name} POST_BUILD
            COMMAND ${CMAKE_COMMAND}
                "-DPANCAKE_TARGET_FILE=$<TARGET_FILE:${target_name}>"
                "-DPANCAKE_TARGET_DIR=$<TARGET_FILE_DIR:${target_name}>"
                "-DPANCAKE_MINGW_BIN=${PANCAKE_MINGW_BIN}"
                -P "${CMAKE_SOURCE_DIR}/cmake/deploy_mingw_runtime.cmake"
            COMMENT "Copying MinGW runtime DLLs next to ${target_name}"
            VERBATIM)
    endif()
endfunction()

set(PANCAKE_HEADERS
    src/communicator.hpp
    src/config.hpp
    src/energy.hpp
    src/equation_of_state.hpp
    src/field.hpp
    src/film_thickness.hpp
    src/fvm.hpp
    src/io.hpp
    src/journal_motion.hpp
    src/linear_system.hpp
    src/mesh.hpp
    src/reynolds.hpp
    src/utils.hpp)

if(PANCAKE_BUILD_SOLVER)
    # Dependencies
    find_package(MPI REQUIRED)
    if(NOT MPIEXEC_EXECUTABLE)
        find_program(MPIEXEC_EXECUTABLE
            NAMES mpiexec.exe mpiexec
            HINTS
                "$ENV{MSMPI_BIN}"
                "C:/Program Files/Microsoft MPI/Bin"
                "C:/Program Files (x86)/Microsoft SDKs/MPI"
            PATH_SUFFIXES Bin bin)
    endif()
    if(MPIEXEC_EXECUTABLE AND NOT MPIEXEC_NUMPROC_FLAG)
        set(MPIEXEC_NUMPROC_FLAG -n)
    endif()

    find_package(PkgConfig QUIET)
    set(PETSC_TARGET "")
    set(PANCAKE_PETSC_PKG_CONFIG_NAMES
        petsc-dmo
        petsc-smo
        petsc-zmo
        petsc-cmo
        PETSc
        petsc)

    if(PkgConfig_FOUND)
        foreach(petsc_pkg_config_name IN LISTS PANCAKE_PETSC_PKG_CONFIG_NAMES)
            if(NOT PETSC_TARGET)
                pkg_check_modules(PETSC QUIET IMPORTED_TARGET ${petsc_pkg_config_name})
                if(PETSC_FOUND)
                    set(PETSC_TARGET PkgConfig::PETSC)
                    message(STATUS "PETSc pkg-config module: ${petsc_pkg_config_name}")
                endif()
            endif()
        endforeach()
    endif()

    if(NOT PETSC_TARGET)
        set(PETSC_DIR "/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real"
            CACHE PATH "PETSc installation directory")
        find_path(PETSC_INCLUDE_DIR petsc.h
            HINTS
                "${PETSC_DIR}/include"
                "${PETSC_DIR}/include/petsc"
                "${PETSC_DIR}/include/petsc/dmo"
                "${PETSC_DIR}/include/petsc/smo"
                "${PETSC_DIR}/include/petsc/zmo"
                "${PETSC_DIR}/include/petsc/cmo"
            REQUIRED)
        find_library(PETSC_LIBRARY
            NAMES petsc_real petsc petsc-dmo petsc-smo petsc-zmo petsc-cmo
            HINTS "${PETSC_DIR}/lib"
            REQUIRED)
        message(STATUS "PETSc include: ${PETSC_INCLUDE_DIR}")
        message(STATUS "PETSc library: ${PETSC_LIBRARY}")
    endif()
    message(STATUS "MPI C++ Compiler: ${MPI_CXX_COMPILER}")

    set(PANCAKE_CORE_SOURCES
        src/communicator.cpp
        src/energy.cpp
        src/field.cpp
        src/film_thickness.cpp
        src/fvm.cpp
        src/io.cpp
        src/journal_motion.cpp
        src/linear_system.cpp
        src/mesh.cpp
        src/reynolds.cpp)

    add_library(pancake_core ${PANCAKE_CORE_SOURCES} ${PANCAKE_HEADERS})
    target_compile_features(pancake_core PUBLIC cxx_std_20)
    target_include_directories(pancake_core PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(pancake_core PUBLIC MPI::MPI_CXX)
    if(MINGW)
        target_compile_definitions(pancake_core PUBLIC
            PETSC_FUNCTION_NAME_CXX=__func__
            PETSC_CXX_RESTRICT=__restrict__)
    endif()
    if(PETSC_TARGET)
        target_link_libraries(pancake_core PUBLIC ${PETSC_TARGET})
    else()
        target_include_directories(pancake_core PUBLIC ${PETSC_INCLUDE_DIR})
        target_link_libraries(pancake_core PUBLIC ${PETSC_LIBRARY})
    endif()

    if(MSVC)
        target_compile_options(pancake_core PRIVATE /W4)
    else()
        target_compile_options(pancake_core PRIVATE -Wall -Wextra)
    endif()

    add_executable(pancake src/main.cpp)
    target_link_libraries(pancake PRIVATE pancake_core)
    if(WIN32)
        enable_language(RC)
        target_sources(pancake PRIVATE src/pancake.rc)
        set_source_files_properties(src/pancake.rc PROPERTIES
            INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/src
            OBJECT_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/pancake.ico")
    endif()
    pancake_copy_default_config(pancake)
    pancake_deploy_mingw_runtime(pancake)

    if(PANCAKE_BUILD_TESTS)
        enable_testing()
        add_subdirectory(tests)
    endif()
else()
    message(STATUS "PANCAKE_BUILD_SOLVER is OFF; skipping MPI/PETSc solver targets.")
endif()

if(WIN32 AND PANCAKE_BUILD_GUI)
    # The ComCtl32 v6 dependency (modern visual styles) is embedded through a
    # manifest resource. MSVC also honours the #pragma in gui_win32.cpp, but
    # MinGW needs the compiled .rc, so wire it in for every Windows compiler.
    enable_language(RC)
    add_executable(pancake_gui WIN32
        src/gui_win32.cpp
        src/pancake_gui.rc
        src/config.hpp)
    target_compile_features(pancake_gui PRIVATE cxx_std_20)
    target_include_directories(pancake_gui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    # windres resolves the manifest path relative to its include dirs.
    set_source_files_properties(src/pancake_gui.rc PROPERTIES
        INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/src
        OBJECT_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/pancake.ico;${CMAKE_CURRENT_SOURCE_DIR}/src/pancake_gui.manifest")
    target_compile_definitions(pancake_gui PRIVATE UNICODE _UNICODE NOMINMAX)
    target_link_libraries(pancake_gui PRIVATE shell32 comctl32 gdi32 comdlg32 ole32)
    # Ship the GUI as a single self-contained pancake_gui.exe: it depends only on
    # Windows system DLLs, so we statically link the C/C++ runtime instead of
    # deploying MinGW runtime DLLs beside it.
    if(MINGW)
        target_link_options(pancake_gui PRIVATE -municode -static)
    elseif(MSVC)
        set_property(TARGET pancake_gui PROPERTY
            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()
    pancake_copy_default_config(pancake_gui)
elseif(PANCAKE_BUILD_GUI)
    message(STATUS "PANCAKE_BUILD_GUI is ON, but the native GUI is Windows-only.")
endif()
