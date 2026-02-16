#include "logger.hpp"
#include <iostream>

namespace Logger {

    // ANSI Escape Codes
    const std::string ANSI_RESET        = "\e[0;0m";
    const std::string ANSI_RED          = "\e[0;91m";
    const std::string ANSI_GREEN        = "\e[0;92m";
    const std::string ANSI_YELLOW       = "\e[0;93m";
    const std::string ANSI_BLUE         = "\e[0;94m";
    const std::string ANSI_CYAN         = "\e[0;96m";
    const std::string ANSI_MAGENTA      = "\e[0;95m";
    const std::string ANSI_BOLD_RED     = "\e[1;91m";
    const std::string ANSI_BOLD_GREEN   = "\e[1;92m";
    const std::string ANSI_BOLD_YELLOW  = "\e[1;93m";
    const std::string ANSI_BOLD_BLUE    = "\e[1;94m";
    const std::string ANSI_BOLD_CYAN    = "\e[1;96m";
    const std::string ANSI_BOLD_MAGENTA = "\e[1;95m";

    // Helper to map enum to string
    std::string get_color_code(Color c) {
        switch (c) {
            case Color::RED:            return ANSI_RED;
            case Color::GREEN:          return ANSI_GREEN;
            case Color::YELLOW:         return ANSI_YELLOW;
            case Color::BLUE:           return ANSI_BLUE;
            case Color::CYAN:           return ANSI_CYAN;
            case Color::MAGENTA:        return ANSI_MAGENTA;
            case Color::BOLD_RED:       return ANSI_BOLD_RED;
            case Color::BOLD_GREEN:     return ANSI_BOLD_GREEN;
            case Color::BOLD_YELLOW:    return ANSI_BOLD_YELLOW;
            case Color::BOLD_BLUE:      return ANSI_BOLD_BLUE;
            case Color::BOLD_CYAN:      return ANSI_BOLD_CYAN;
            case Color::BOLD_MAGENTA:   return ANSI_BOLD_MAGENTA;
            default:                return ANSI_RESET;
        }
    }

    void print(const std::string& msg, Color color, int indent_level, bool rank_zero_only) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // Filter: If rank_zero_only is true, only Rank 0 prints.
        if (rank_zero_only && rank != 0) {
            return;
        }

        // Build prefix
        std::string prefix = "";
        for (int i = 0; i < indent_level; ++i) {
            prefix += "  "; // 2 spaces per indent
        }

        // Print to stdout
        // Note: std::endl flushes the buffer, which is important for MPI logs
        std::cout << get_color_code(color) << prefix << msg << ANSI_RESET << std::endl;
    }
}
