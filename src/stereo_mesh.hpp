#pragma once

#include <cmath>
#include <mpi.h>
#include <vector>

#include "bubble_config.hpp"

/// Minimal 3D vector with arithmetic operations.
struct Vec3 {
    double x, y, z;
    Vec3 operator+(const Vec3& b) const { return {x+b.x, y+b.y, z+b.z}; }
    Vec3 operator-(const Vec3& b) const { return {x-b.x, y-b.y, z-b.z}; }
    Vec3 operator*(double s)      const { return {x*s,   y*s,   z*s};   }
    double dot(const Vec3& b)     const { return x*b.x + y*b.y + z*b.z; }
    Vec3   cross(const Vec3& b)   const {
        return {y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x};
    }
    double norm()      const { return std::sqrt(x*x + y*y + z*z); }
    Vec3   normalized() const { double n = norm(); return {x/n, y/n, z/n}; }
};

/// Structured 2D grid on a hemispherical surface via inverse stereographic projection.
///
/// The 2D domain is the square [-L_box, L_box]^2 with n_u x n_v uniform cells.
/// Each cell center (u_i, v_j) maps to a 3D point on the hemisphere via:
///   x = 2R^2 u / (r^2 + R^2),   y = 2R^2 v / (r^2 + R^2),
///   z = R (R^2 - r^2) / (r^2 + R^2),   r^2 = u^2 + v^2.
/// Cells with r > R_bubble map below the equator and are handled by the Mask.
///
/// Domain decomposition: 1D in u (rows), matching the bearing pattern.
/// All v-indices are global (v is not decomposed).
///
/// Geometry arrays are laid out as [i * n_v_global + j] for
/// local i in [0, n_u_local), j in [0, n_v_global).
class StereoMesh {
public:
    int rank, size;
    int n_u_global, n_v_global;
    int n_u_local, n_v_local;   // n_v_local == n_v_global (v not decomposed)
    int offset_u;
    double R, L_box, du, dv;

    // Per-cell geometry (physical cells only, indexed [i * n_v_global + j])
    std::vector<double> cell_area;   // 3D projected cell area [m^2]
    std::vector<double> face_len_e;  // East face 3D chord length (v-aligned) [m]
    std::vector<double> face_len_n;  // North face 3D chord length (u-aligned) [m]
    std::vector<double> d_u;         // Center-to-center 3D distance (i,j)->(i+1,j) [m]
    std::vector<double> d_v;         // Center-to-center 3D distance (i,j)->(i,j+1) [m]
    std::vector<double> cell_nx, cell_ny, cell_nz;  // Outward unit normal at cell centre
    std::vector<double> g_u, g_v;    // Surface gravity in local (u,v) tangent directions [m/s^2]

    StereoMesh(const BubbleConfig& cfg, MPI_Comm comm = MPI_COMM_WORLD);

    /// Inverse stereographic projection: (u, v) -> 3D point on hemisphere.
    Vec3 project(double u, double v) const;

    /// 2D coordinate of local cell centre in u direction.
    double cell_u(int local_i) const { return -L_box + (offset_u + local_i + 0.5) * du; }
    /// 2D coordinate of cell centre in v direction.
    double cell_v(int j)       const { return -L_box + (j + 0.5) * dv; }

    /// 2D coordinate of the left vertex of local cell i in u (i.e., face at i - 1/2).
    double vertex_u(int local_i) const { return -L_box + (offset_u + local_i) * du; }
    /// 2D coordinate of the bottom vertex of cell j in v.
    double vertex_v(int j)       const { return -L_box + j * dv; }

    int global_u(int local_i) const { return offset_u + local_i; }
    int idx(int i, int j)     const { return i * n_v_global + j; }

    /// East face chord length for any global u index (computes analytically if not local).
    double east_face_len(int global_u, int j) const;
    /// 3D centre-to-centre distance from (global_u, j) to (global_u+1, j).
    double d_u_east(int global_u, int j) const;

private:
    void compute_geometry(double gx, double gy, double gz);
};
