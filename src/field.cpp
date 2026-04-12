#include "field.hpp"

Field::Field(std::string name, const Mesh& mesh, int ghost_layers, GridLocation loc)
    : name(name), location(loc), n_ghost(ghost_layers), mesh_ref(mesh) {
    n_theta_phys = mesh.n_theta_local;
    n_z_phys = mesh.n_z_local;
    
    if (location == GridLocation::FACE_THETA) {
        // One more face in theta direction
        // For MPI, the extra face is at the end of the global domain.
        // But since it's periodic, face N is same as face 0.
        // However, for local indexing within a rank, it's often easier to 
        // have n+1 faces if we want to avoid ghost lookups for the right face.
        // Actually, if it's periodic, n faces are enough.
        // But the user asked for staggered grid, and standard implementation 
        // usually has N faces for periodic N cells.
        // Let's stick to N for now if it's simpler, or N+1 if we want to be explicit.
        // Given Communicator logic, N faces works perfectly for periodic theta.
    } else if (location == GridLocation::FACE_Z) {
        // One more face in z direction (non-periodic)
        n_z_phys = mesh.n_z_local + 1;
    }

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
