# Source File Audit

Status of each source file in the bubble solver.

| File | Role | Notes |
|------|------|-------|
| `bubble_config.hpp/cpp` | All solver parameters + `ConvectionScheme` | Loaded from `bubble_config.txt` at runtime |
| `stereo_mesh.hpp/cpp` | Inverse stereographic grid + precomputed 3D geometry | 1D MPI decomposition in u |
| `mask.hpp/cpp` | Circular domain mask + rim BC helpers | `active_at()` helper works for any global index |
| `communicator.hpp/cpp` | MPI ghost exchange, open (non-periodic) | Zero-gradient fill at all domain edges |
| `bubble_io.hpp/cpp` | VTK unstructured output (vtu/pvtu/pvd) | One `.vtu` per rank; rank 0 writes `.pvtu` and `.pvd` |
| `field.hpp/cpp` | Scalar fields with ghost layers | `operator()(i,j)` handles ghost offset arithmetic |
| `bubble_linear_system.hpp/cpp` | PETSc 5-point sparse system | `apply_mask()` does row replacement for masked/rim cells |
| `surface_operators.hpp/cpp` | Surface FVM: ddt, laplacian, divergence, add_source, gradient | Uses 3D geometry from StereoMesh throughout |
| `thickness_transport.hpp/cpp` | Gravity velocity, face fluxes, h transport solve | Phase 3 complete |
| `utils.hpp` | Rank-0 coloured logger | Header-only |
| `bubble_main.cpp` | Entry point and time loop | Phases 1–3 wired; Phases 4–6 have TODO stubs |
