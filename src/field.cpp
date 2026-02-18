#include "field.hpp"

Field::Field(std::string name, const Mesh& mesh, int ghost_layers, GridLocation loc)
    : name(name), location(loc), n_ghost(ghost_layers), mesh_ref(mesh) {
    n_theta_phys = mesh.n_theta_local;
    n_z_phys = mesh.n_z_local;
    stride_theta = n_theta_phys + 2 * n_ghost;
    data.resize(stride_theta * (n_z_phys + 2 * n_ghost), 0.0);
    old_data = data;
}

void Field::fill(double value) { std::fill(data.begin(), data.end(), value); }

Field& Fields::add(const std::string& name, const Mesh& mesh, int n_ghost, GridLocation loc) {
    if(storage.find(name) != storage.end()) throw std::runtime_error("Field exists: " + name);
    storage[name] = std::make_unique<Field>(name, mesh, n_ghost, loc);
    return *storage[name];
}

Field& Fields::operator[](const std::string& name) { return *storage.at(name); }
const Field& Fields::operator[](const std::string& name) const { return *storage.at(name); }
