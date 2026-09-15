#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

enum class GridLocation { CENTER, FACE_THETA, FACE_Z };

/// Cell-centered (or face-staggered) scalar field with ghost layers.
/// Indexing: operator()(i, j) with i in [-n_ghost, n_u_phys + n_ghost)
///                                   j in [-n_ghost, n_v_phys + n_ghost).
/// Physical cells occupy indices [0, n_u_phys) x [0, n_v_phys).
/// In bearing context: u = theta, v = z.
/// In bubble context:  u = stereographic u, v = stereographic v.
class Field {
public:
    std::string name;
    GridLocation location;
    int n_theta_phys, n_z_phys, n_ghost, stride_theta;
    std::vector<double> data;
    std::vector<double> old_data;

    /// Primary constructor: takes raw physical dimensions.
    Field(std::string name, int n_u, int n_v, int n_ghost = 2,
          GridLocation loc = GridLocation::CENTER);

    inline double& operator()(int i, int j) {
        return data[(j + n_ghost) * stride_theta + (i + n_ghost)];
    }
    inline double operator()(int i, int j) const {
        return data[(j + n_ghost) * stride_theta + (i + n_ghost)];
    }
    inline double old(int i, int j) const {
        return old_data[(j + n_ghost) * stride_theta + (i + n_ghost)];
    }

    void fill(double value);
    void store_old_time() { old_data = data; }
};

/// Named container for Fields. Provides safe lookup and construction.
class Fields {
private:
    std::map<std::string, std::unique_ptr<Field>> storage;
public:
    /// Add a field with explicit physical dimensions.
    Field& add(const std::string& name, int n_u, int n_v,
               int n_ghost = 2, GridLocation loc = GridLocation::CENTER);

    Field& operator[](const std::string& name);
    const Field& operator[](const std::string& name) const;
    bool has(const std::string& name) const { return storage.count(name) > 0; }

    auto begin() { return storage.begin(); }
    auto end()   { return storage.end(); }
    auto begin() const { return storage.cbegin(); }
    auto end()   const { return storage.cend(); }
};
