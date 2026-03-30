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

## Test Plan

All tests must run on a workstation without cluster access (use a coarse mesh).
C++ tests are registered with CTest; Python tests use `unittest`.
Run the full suite with `ctest --output-on-failure` from the build directory,
and `python -m pytest tests/` for Python tests.

### Phase 1 — Mesh (`tests/test_mesh.py`)
| Check | Pass criterion |
|-------|---------------|
| `generate_mesh.py` runs without error | exit code 0 |
| Output `.msh` file exists | file present |
| Node count in expected range | 3 000 – 25 000 nodes |
| Triangle count in expected range | 5 000 – 50 000 elements |
| All 4 physical boundary tags present | tags 1–4 found in mesh |
| No degenerate elements (quality check via Gmsh API) | min quality > 0.1 |

### Phase 2 — FEM system initialisation (`tests/test_spaces.cpp`)
| Check | Pass criterion |
|-------|---------------|
| libMesh initialises without error | no exception |
| Mesh loads from `.msh` file | `mesh.n_elem() > 0` |
| P2 velocity space has more DOFs than P1 pressure space | `n_u_dofs > n_p_dofs` |
| Inlet Dirichlet DOFs are marked | at least 1 constrained DOF on boundary 1 |
| No-slip DOFs on cylinder are marked | at least 1 constrained DOF on boundary 4 |

### Phase 3 — Stokes initialisation solve (`tests/test_stokes.cpp`)
| Check | Pass criterion |
|-------|---------------|
| Stokes linear solve converges | KSP converged reason > 0 |
| KSP iteration count | < 300 |
| Linear residual after solve | < 1e-8 |
| Velocity is zero at no-slip boundaries | max |u| on walls/cylinder < 1e-10 |
| Pressure has zero mean (or pinned reference) | solver does not diverge |

### Phase 4 — Output and drag/lift (`tests/test_output.cpp`)
| Check | Pass criterion |
|-------|---------------|
| ExodusII file is written | file `results/test.e` exists |
| File can be re-read by libMesh ExodusII reader | no exception on read-back |
| Drag force on cylinder is finite and non-NaN | `std::isfinite(F_D)` |
| Lift force is near zero for symmetric Stokes flow | `|C_L| < 1e-6` |

### Phase 5 — CMake build (`CTest` infrastructure)
| Check | Pass criterion |
|-------|---------------|
| `cmake` configure succeeds | exit code 0 |
| `cmake --build` compiles all targets | exit code 0, binary present |
| All CTest tests discovered | `ctest -N` lists ≥ 3 tests |
| Solver binary runs with `--help` | exit code 0 or 1, no crash |

### Phase 6 — Python driver (`tests/test_driver.py`)
| Check | Pass criterion |
|-------|---------------|
| CLI argument parsing works | all args parsed without error |
| SLURM script is generated with correct structure | `#SBATCH` directives present |
| Mesh generation is invoked correctly | `generate_mesh.py` called with right args |
| Dry-run mode does not submit or execute | no side effects |

### Phase 7 — Validation (`tests/test_validation.py`)
| Check | Pass criterion |
|-------|---------------|
| Re = 20 steady solve converges | Newton converged, ≤ 20 iterations |
| Drag coefficient (Re = 20) | C_D ∈ [5.57, 5.59] (Schafer-Turek interval) |
| Lift coefficient (Re = 20) | \|C_L\| < 0.02 |
| Lift coefficient (Re = 5) | \|C_L\| < 1e-4 (symmetric flow) |
| ExodusII output is written | file present and non-empty |

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
│   ├── benchmark_targets.md         # Expected C_D, C_L values
│   └── conversation_log.md          # Session-by-session decision record
├── src/
│   ├── CMakeLists.txt               # Build system (includes tests/)
│   ├── params.h                     # Centralized numerical parameters
│   ├── channel_flow_system.h        # FEMSystem subclass declaration
│   ├── channel_flow_system.cpp      # Weak form assembly + BCs
│   └── main.cpp                     # Entry point
├── tests/
│   ├── CMakeLists.txt               # CTest registration for C++ tests
│   ├── test_mesh.py                 # Phase 1: mesh validity checks
│   ├── test_spaces.cpp              # Phase 2: DOF space initialisation
│   ├── test_stokes.cpp              # Phase 3: Stokes linear solve
│   ├── test_output.cpp              # Phase 4: ExodusII I/O + drag/lift
│   ├── test_driver.py               # Phase 6: Python driver logic
│   └── test_validation.py           # Phase 7: benchmark C_D, C_L
├── scripts/
│   └── run_simulation.py            # Python driver
└── results/                         # (git-ignored) simulation output
```
