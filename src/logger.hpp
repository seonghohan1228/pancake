#pragma once
#include <string>
#include <mpi.h>

namespace Logger {

    enum class Color {
        RESET,
        RED,
        GREEN,
        YELLOW,
        BLUE,
        CYAN,
        MAGENTA,
        BOLD_RED,
        BOLD_GREEN,
        BOLD_YELLOW,
        BOLD_BLUE,
        BOLD_CYAN,
        BOLD_MAGENTA,
    };

    // Prints a formatted message to stdout with MPI rank filtering.
    void print(const std::string& msg,
               Color color = Color::RESET,
               int indent_level = 0,
               bool rank_zero_only = true);
}
