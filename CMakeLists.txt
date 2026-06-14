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

set(PANCAKE_CASE_CONFIGS
    config_r290_pz68.txt
    config_r290_pz68_quick.txt)

foreach(PANCAKE_CASE_CONFIG IN LISTS PANCAKE_CASE_CONFIGS)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/${PANCAKE_CASE_CONFIG}
        ${CMAKE_CURRENT_BINARY_DIR}/${PANCAKE_CASE_CONFIG}
        COPYONLY)
endforeach()

function(pancake_copy_default_config target_name)
    add_custom_command(TARGET ${target_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/config.txt"
            "$<TARGET_FILE_DIR:${target_name}>/config.txt"
        COMMENT "Copying config.txt next to target executable"
        VERBATIM)
endfunction()

function(pancake_copy_case_configs target_name)
    foreach(PANCAKE_CASE_CONFIG IN LISTS PANCAKE_CASE_CONFIGS)
        add_custom_command(TARGET ${target_name} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_CURRENT_SOURCE_DIR}/${PANCAKE_CASE_CONFIG}"
                "$<TARGET_FILE_DIR:${target_name}>/${PANCAKE_CASE_CONFIG}"
            COMMENT "Copying ${PANCAKE_CASE_CONFIG} next to target executable"
            VERBATIM)
    endforeach()
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
    src/diagnostics.hpp
    src/energy.hpp
    src/equation_of_state.hpp
    src/field.hpp
    src/fluid_properties.hpp
    src/film_thickness.hpp
    src/fvm.hpp
    src/io.hpp
    src/journal_motion.hpp
    src/linear_system.hpp
    src/mesh.hpp
    src/property_table_io.hpp
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
        src/diagnostics.cpp
        src/energy.cpp
        src/field.cpp
        src/fluid_properties.cpp
        src/film_thickness.cpp
        src/fvm.cpp
        src/gas_transport.cpp
        src/io.cpp
        src/journal_motion.cpp
        src/linear_system.cpp
        src/mesh.cpp
        src/property_table_io.cpp
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
    pancake_copy_case_configs(pancake)
    pancake_deploy_mingw_runtime(pancake)

    if(PANCAKE_BUILD_TESTS)
        enable_testing()
        add_subdirectory(tests)
    endif()
else()
    message(STATUS "PANCAKE_BUILD_SOLVER is OFF; skipping MPI/PETSc solver targets.")
endif()

if(WIN32 AND PANCAKE_BUILD_GUI)
    enable_language(RC)
    include(FetchContent)

    # Dear ImGui (docking branch) + ImPlot, pinned releases fetched at
    # configure time. No system packages, no DLLs - everything is compiled
    # into pancake_gui.exe.
    FetchContent_Declare(imgui
        URL https://github.com/ocornut/imgui/archive/refs/tags/v1.91.9b-docking.zip)
    FetchContent_Declare(implot
        URL https://github.com/epezent/implot/archive/refs/tags/v0.17.zip)
    FetchContent_MakeAvailable(imgui implot)

    add_library(pancake_imgui STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_win32.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_dx11.cpp
        ${implot_SOURCE_DIR}/implot.cpp
        ${implot_SOURCE_DIR}/implot_items.cpp)
    target_compile_features(pancake_imgui PUBLIC cxx_std_20)
    target_include_directories(pancake_imgui PUBLIC
        ${imgui_SOURCE_DIR}
        ${imgui_SOURCE_DIR}/backends
        ${implot_SOURCE_DIR})
    target_compile_definitions(pancake_imgui PUBLIC
        IMGUI_DISABLE_OBSOLETE_FUNCTIONS
        IMGUI_DISABLE_DEMO_WINDOWS
        IMGUI_IMPL_WIN32_DISABLE_GAMEPAD
        UNICODE _UNICODE NOMINMAX)

    # Embedded assets: fonts are converted to C arrays at build time, so the
    # exe reads zero asset files from disk at startup.
    set(PANCAKE_GUI_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
    file(MAKE_DIRECTORY ${PANCAKE_GUI_GENERATED_DIR})
    function(pancake_embed_binary input_file symbol header_name)
        add_custom_command(
            OUTPUT ${PANCAKE_GUI_GENERATED_DIR}/${header_name}
            COMMAND ${CMAKE_COMMAND}
                -DEMBED_INPUT=${input_file}
                -DEMBED_OUTPUT=${PANCAKE_GUI_GENERATED_DIR}/${header_name}
                -DEMBED_SYMBOL=${symbol}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_binary.cmake
            DEPENDS ${input_file} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_binary.cmake
            COMMENT "Embedding ${header_name}"
            VERBATIM)
    endfunction()
    pancake_embed_binary(${CMAKE_CURRENT_SOURCE_DIR}/third_party/fonts/Roboto-Regular.ttf
        g_font_roboto_regular embedded_roboto_regular.h)
    pancake_embed_binary(${CMAKE_CURRENT_SOURCE_DIR}/third_party/fonts/Roboto-Medium.ttf
        g_font_roboto_medium embedded_roboto_medium.h)
    pancake_embed_binary(${CMAKE_CURRENT_SOURCE_DIR}/third_party/fonts/JetBrainsMono-Regular.ttf
        g_font_jetbrains_mono embedded_jetbrains_mono.h)

    add_executable(pancake_gui WIN32
        src/gui/main_gui.cpp
        src/gui/gui_app.cpp
        src/gui/panels_setup.cpp
        src/gui/panels_results.cpp
        src/gui/solver_runner.cpp
        src/gui/results.cpp
        src/gui/schema.cpp
        src/gui/exports.cpp
        src/gui/settings.cpp
        src/gui/theme.cpp
        src/gui/pancake_gui.rc
        ${PANCAKE_GUI_GENERATED_DIR}/embedded_roboto_regular.h
        ${PANCAKE_GUI_GENERATED_DIR}/embedded_roboto_medium.h
        ${PANCAKE_GUI_GENERATED_DIR}/embedded_jetbrains_mono.h)
    target_compile_features(pancake_gui PRIVATE cxx_std_20)
    target_include_directories(pancake_gui PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/gui
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/stb
        ${PANCAKE_GUI_GENERATED_DIR})
    # windres resolves resource paths relative to its include dirs.
    set_source_files_properties(src/gui/pancake_gui.rc PROPERTIES
        INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/src
        OBJECT_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/pancake.ico;${CMAKE_CURRENT_SOURCE_DIR}/src/gui/pancake_gui.manifest")
    target_link_libraries(pancake_gui PRIVATE pancake_imgui
        d3d11 dxgi d3dcompiler dwmapi user32 gdi32 shell32 ole32 comdlg32)

    # Warnings-as-errors for first-party GUI sources only (not third-party).
    if(MSVC)
        target_compile_options(pancake_gui PRIVATE /W4 /WX)
    else()
        target_compile_options(pancake_gui PRIVATE -Wall -Wextra -Werror)
    endif()

    # Single statically linked exe: only Windows-shipped DLLs in the import
    # table (no vcruntime/msvcp, no MinGW runtime DLLs).
    if(MINGW)
        target_link_options(pancake_gui PRIVATE -municode -static)
    elseif(MSVC)
        set_property(TARGET pancake_gui PROPERTY
            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
        set_property(TARGET pancake_imgui PROPERTY
            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()
    # Convenience for local runs; the exe itself does not require these files.
    pancake_copy_default_config(pancake_gui)
    pancake_copy_case_configs(pancake_gui)
elseif(PANCAKE_BUILD_GUI)
    message(STATUS "PANCAKE_BUILD_GUI is ON, but the native GUI is Windows-only.")
endif()
