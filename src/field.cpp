#include "field.hpp"

Field::Field(std::string name, int n_u, int n_v, int n_ghost, GridLocation loc)
    : name(name), location(loc), n_ghost(n_ghost)
{
    n_theta_phys = n_u;
    n_z_phys     = n_v;

    if (loc == GridLocation::FACE_Z) {
        n_z_phys = n_v + 1;  // One extra face in v direction (non-periodic)
    }

    stride_theta = n_theta_phys + 2 * n_ghost;
    data.resize(stride_theta * (n_z_phys + 2 * n_ghost), 0.0);
    old_data = data;
}

void Field::fill(double value) { std::fill(data.begin(), data.end(), value); }

Field& Fields::add(const std::string& name, int n_u, int n_v, int n_ghost, GridLocation loc) {
    if (storage.count(name)) throw std::runtime_error("Field already exists: " + name);
    storage[name] = std::make_unique<Field>(name, n_u, n_v, n_ghost, loc);
    return *storage[name];
}

Field& Fields::operator[](const std::string& name) { return *storage.at(name); }
const Field& Fields::operator[](const std::string& name) const { return *storage.at(name); }
