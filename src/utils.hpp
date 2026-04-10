#pragma once

#include <iostream>
#include <mpi.h>
#include <string>

#include "field.hpp"

namespace Utils {
    // Logger (console output)
    enum class Color { RESET, RED, GREEN, YELLOW, BLUE, CYAN, MAGENTA };

    inline std::string get_code(Color c) {
        switch(c) {
            case Color::RED: return "\033[31m";
            case Color::GREEN: return "\033[32m";
            case Color::YELLOW: return "\033[33m";
            case Color::BLUE: return "\033[34m";
            case Color::CYAN: return "\033[36m";
            default: return "\033[0m";
        }
    }

    inline void log(const std::string& msg, Color color = Color::RESET, int indent = 0) {
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            std::string prefix(indent * 2, ' ');
            std::cout << get_code(color) << prefix << msg << get_code(Color::RESET) << std::endl;
        }
    }

    // --- Boundary Conditions ---
    enum class EndCondition { FIXED, ZERO_GRADIENT };

    inline void apply_z_boundaries(Field& f, EndCondition type, double value = 0.0) {
        int n_theta = f.n_theta_phys;
        int n_z = f.n_z_phys;

        for (int i = 0; i < n_theta; ++i) {
            // Bottom Boundary (z=0, ghost index -1)
            double val_0 = f(i, 0);
            f(i, -1) = (type == EndCondition::FIXED) ? value : val_0;

            // Top Boundary (z=L, ghost index nz)
            double val_last = f(i, n_z - 1);
            f(i, n_z) = (type == EndCondition::FIXED) ? value : val_last;
        }
    }
}
