# Implementation Plan — 2D Channel Flow Past a Cylinder

## Architecture
- **Solver**: C++ with libMesh (FEM) + PETSc (linear algebra/nonlinear solver)
- **Driver**: Python (`scripts/run_simulation.py`) — mesh generation, job submission, post-processing
- **Benchmark**: Schafer & Turek (1996) DFG 2D-1

## Phases

| # | Phase | Status |
|---|-------|--------|
| 0 | Project scaffolding (dirs, PLAN.md, notes, .gitignore) | ✅ Done |
| 1 | Gmsh mesh generation | ✅ Done |
| 2 | C++ solver — libMesh FEMSystem (Taylor-Hood P2/P1) | ⬜ Todo |
| 3 | Solver configuration — PETSc SNES + FGMRES + BoomerAMG | ⬜ Todo |
| 4 | Output — ExodusII + drag/lift post-processing | ⬜ Todo |
| 5 | CMake build system | ⬜ Todo |
| 6 | Python driver | ⬜ Todo |
| 7 | Validation against benchmark | ⬜ Todo |

## Key Decisions

| Decision | Choice | Justification |
|----------|--------|---------------|
| FEM framework | libMesh | Pure C++ solver; PETSc + Trilinos backends available |
| Elements | Taylor-Hood P2/P1 | Inf-sup stable (Brezzi 1974); standard for Navier-Stokes |
| Weak form | Standard Galerkin | No stabilization needed at Re ≈ 5 |
| Advection form | Skew-symmetric | Energy stable |
| Nonlinear solver | PETSc SNES Newton w/ backtracking line search | Quadratic convergence near solution |
| Linear solver | FGMRES (restart=100) | Flexible Krylov for nonsymmetric saddle-point system |
| Preconditioner | BoomerAMG (Hypre) | Optimal complexity for elliptic-dominated systems at low Re |
| Initialization | Stokes → Navier-Stokes | Avoids Newton divergence from zero initial guess |
| Output format | ExodusII (.e) | Native to libMesh; ParaView-compatible |
| Mesh | Gmsh unstructured triangles | h=0.05 global, h=0.01 near cylinder |

## Domain Parameters (Schafer-Turek 2D-1)

```
Channel:  [0, 2.2] × [0, 0.41]   (meters)
Cylinder: center (0.2, 0.2), radius 0.05
U_max:    0.3 m/s  (mean inlet velocity U_mean = 0.2 m/s)
nu:       1e-3 m²/s  (kinematic viscosity, water-like)
Re:       U_mean * D / nu = 0.2 * 0.1 / 1e-3 = 20  (Schafer-Turek)
          For Re=5: U_mean = 0.05 m/s, U_max = 0.075 m/s
```

## File Map

```
.
├── CLAUDE.md                        # Project guidelines
├── PLAN.md                          # This file — living tracker
├── .gitignore
├── meshes/
│   ├── channel.geo                  # Gmsh geometry (parametric)
│   └── generate_mesh.py             # Mesh generation script
├── notes/
│   ├── solver_parameters.md         # Literature references for numerical choices
│   └── benchmark_targets.md         # Expected C_D, C_L values
├── src/
│   ├── CMakeLists.txt               # Build system
│   ├── channel_flow_system.h        # FEMSystem subclass declaration
│   ├── channel_flow_system.cpp      # Weak form assembly + BCs
│   └── main.cpp                     # Entry point
├── scripts/
│   └── run_simulation.py            # Python driver
└── results/                         # (git-ignored) simulation output
```
