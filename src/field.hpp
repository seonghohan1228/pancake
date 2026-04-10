#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "mesh.hpp"

enum class GridLocation { CENTER, FACE_THETA, FACE_Z };

class Field {
public:
    std::string name;
    GridLocation location;
    int n_theta_phys, n_z_phys, n_ghost, stride_theta;
    std::vector<double> data;
    std::vector<double> old_data; // For time derivatives
    const Mesh& mesh_ref;

    Field(std::string name, const Mesh& mesh, int ghost_layers = 2, GridLocation loc = GridLocation::CENTER);

    inline double& operator()(int i, int j) { return data[(j + n_ghost) * stride_theta + (i + n_ghost)]; }      // Read/write
    inline double operator()(int i, int j) const { return data[(j + n_ghost) * stride_theta + (i + n_ghost)]; } // Read-only

    void fill(double value);
    void store_old_time() { old_data = data; }
};

class Fields {
private:
    std::map<std::string, std::unique_ptr<Field>> storage;
public:
    Field& add(const std::string& name, const Mesh& mesh, int ghost_layers = 2, GridLocation loc = GridLocation::CENTER);
    Field& operator[](const std::string& name);
    const Field& operator[](const std::string& name) const;

    // For iterating
    auto begin() { return storage.begin(); }
    auto end() { return storage.end(); }
};
