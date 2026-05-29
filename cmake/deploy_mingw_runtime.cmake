if(NOT DEFINED PANCAKE_TARGET_FILE)
    message(FATAL_ERROR "PANCAKE_TARGET_FILE is required")
endif()

if(NOT DEFINED PANCAKE_TARGET_DIR)
    message(FATAL_ERROR "PANCAKE_TARGET_DIR is required")
endif()

if(NOT DEFINED PANCAKE_MINGW_BIN)
    message(FATAL_ERROR "PANCAKE_MINGW_BIN is required")
endif()

find_program(PANCAKE_LDD_EXECUTABLE
    NAMES ldd
    HINTS
        "${PANCAKE_MINGW_BIN}"
        "${PANCAKE_MINGW_BIN}/../../usr/bin")

if(NOT PANCAKE_LDD_EXECUTABLE)
    message(WARNING "Could not find ldd; MinGW runtime DLLs were not copied.")
    return()
endif()

execute_process(
    COMMAND "${PANCAKE_LDD_EXECUTABLE}" "${PANCAKE_TARGET_FILE}"
    RESULT_VARIABLE ldd_result
    OUTPUT_VARIABLE ldd_output
    ERROR_VARIABLE ldd_error
    OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT ldd_result EQUAL 0)
    message(WARNING "ldd failed for ${PANCAKE_TARGET_FILE}: ${ldd_error}")
    return()
endif()

string(REPLACE "\n" ";" ldd_lines "${ldd_output}")
set(copied_dlls "")

foreach(line IN LISTS ldd_lines)
    if(line MATCHES "=>[ \t]+([^ \t\r\n]+)")
        set(dll_path "${CMAKE_MATCH_1}")
    else()
        continue()
    endif()

    set(source_path "")
    if(dll_path MATCHES "^/mingw64/bin/")
        get_filename_component(dll_name "${dll_path}" NAME)
        set(source_path "${PANCAKE_MINGW_BIN}/${dll_name}")
    elseif(dll_path MATCHES "^[A-Za-z]:/" AND EXISTS "${dll_path}")
        file(TO_CMAKE_PATH "${dll_path}" source_path)
    endif()

    if(NOT source_path)
        continue()
    endif()

    file(TO_CMAKE_PATH "${PANCAKE_MINGW_BIN}" mingw_bin_path)
    if(NOT source_path MATCHES "^${mingw_bin_path}/")
        continue()
    endif()

    get_filename_component(dll_name "${source_path}" NAME)
    list(FIND copied_dlls "${dll_name}" dll_index)
    if(NOT dll_index EQUAL -1)
        continue()
    endif()

    file(COPY_FILE
        "${source_path}"
        "${PANCAKE_TARGET_DIR}/${dll_name}"
        ONLY_IF_DIFFERENT)
    list(APPEND copied_dlls "${dll_name}")
endforeach()

if(copied_dlls)
    list(JOIN copied_dlls ", " copied_dlls_text)
    message(STATUS "Copied MinGW runtime DLLs for ${PANCAKE_TARGET_FILE}: ${copied_dlls_text}")
endif()
