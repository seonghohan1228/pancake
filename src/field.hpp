#pragma once
#include <vector>
#include <string>
#include "mesh.hpp"

enum class GridLocation {
    CENTER,
    FACE_THETA,
    FACE_Z
};

class Field {
public:
    // --- Metadata ---
    std::string name;
    GridLocation location;

    // --- Dimensions ---
    int n_theta_phys;
    int n_z_phys;
    int ng; // Number of ghost layers
    int stride_theta;

    // --- Data Storage ---
    std::vector<double> data;

    // --- Constructor ---
    Field(std::string name, const Mesh& mesh, int ghost_layers = 2, GridLocation loc = GridLocation::CENTER);

    // --- Accessors ---

    // Read/Write
    inline double& operator()(int i, int j) {
        int idx = (j + ng) * stride_theta + (i + ng);
        return data[idx];
    }

    // Read-only
    inline double operator()(int i, int j) const {
        int idx = (j + ng) * stride_theta + (i + ng);
        return data[idx];
    }

    // --- Utilities ---
    void fill(double value);
    double* data_ptr() { return data.data(); }
};

// Operator overloads
// Field + Field
inline Field operator+(const Field& lhs, const Field& rhs) {
    Field result = lhs;
    result.name = "(" + lhs.name + "+" + rhs.name + ")";
    for (size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] + rhs.data[i];
    }
    return result;
}

// Field - Field
inline Field operator-(const Field& lhs, const Field& rhs) {
    Field result = lhs;
    result.name = "(" + lhs.name + "-" + rhs.name + ")";
    for (size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] - rhs.data[i];
    }
    return result;
}

// Field * Field (Element-wise)
inline Field operator*(const Field& lhs, const Field& rhs) {
    Field result = lhs;
    result.name = "(" + lhs.name + "*" + rhs.name + ")";
    for (size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] * rhs.data[i];
    }
    return result;
}

// Field * Scalar
inline Field operator*(const Field& lhs, double scalar) {
    Field result = lhs;
    for (double& val : result.data) val *= scalar;
    return result;
}

// Scalar * Field
inline Field operator*(double scalar, const Field& rhs) {
    return rhs * scalar;
}
