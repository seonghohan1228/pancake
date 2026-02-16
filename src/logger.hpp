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

    /**
     * @brief Prints a formatted message to stdout with MPI rank filtering.
     * * @param msg The message to print.
     * @param color The color of the text.
     * @param indent_level Number of indentations (2 spaces per level).
     * @param arrow If true, adds "-> " before the message (after indent).
     * @param rank_zero_only If true, only the root process prints (default: true).
     */
    void print(const std::string& msg,
               Color color = Color::RESET,
               int indent_level = 0,
               bool rank_zero_only = true);
}
