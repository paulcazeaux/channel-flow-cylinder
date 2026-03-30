# Project Guidelines: 2D Channel Flow Past a Cylinder

## Problem Description
Simulate 2D incompressible viscous flow in a channel past a circular cylinder.
- Reynolds number: Re ≈ 5 (laminar, steady-state eddy regime)
- Expected physics: symmetric recirculation eddies behind the cylinder
- Reference: Schafer & Turek (1996) DFG benchmark "Flow around a cylinder" (2D-1 case)

## Boundary Conditions
- Inlet (left): parabolic velocity profile, u = U_max * 4y(H-y)/H², v = 0
- Outlet (right): stress-free / do-nothing (natural Neumann)
- Top/bottom walls: no-slip (u = v = 0)
- Cylinder surface: no-slip (u = v = 0)

## Discretization & Solver Parameters (from literature)
- Spatial: Taylor-Hood P2/P1 finite elements (or Q2/Q1) — inf-sup stable
- Mesh: unstructured, finer near cylinder; h ~ 0.02–0.05 near cylinder
- Nonlinear solver: Picard iteration or Newton (Newton preferred for Re > 1)
- Linear solver: GMRES with ILU or AMG preconditioner
- Convergence: residual tolerance 1e-6 (nonlinear), 1e-8 (linear)

## Libraries & HPC Environment
- **libMesh** for FEM discretization (pure C++, supports PETSc and Trilinos backends)
- **PETSc** for linear algebra, KSP/PC solvers, and nonlinear solver (SNES)
- **Trilinos** available as alternative/complement (ML/MueLu for AMG, Belos for Krylov)
- Do NOT implement linear solvers, preconditioners, or mesh routines from scratch
- Load libraries via `module load` on the cluster; use placeholders until module names are confirmed:
  ```
  module load gcc/X.Y openmpi/X.Y petsc/X.Y libmesh/X.Y
  ```

## Code Structure
- `src/` — C++ solver source files
- `scripts/` — Python driver (`run_simulation.py`) for setup, job submission, post-processing
- `meshes/` — geometry and mesh files (Gmsh `.geo` / `.msh`)
- `results/` — output (ExodusII `.e` files for ParaView; git-ignored)
- `notes/` — literature justification for numerical choices and benchmark targets

## Code Style
- All functions and classes must have concise Doxygen comments
- Comments only where logic is non-obvious
- **Source files must not exceed 300 lines** — split by logical concern if needed
- **All numerical parameters** (mesh sizes, tolerances, solver flags, physical constants) must
  live in a single centralized header `src/params.h`; no magic numbers elsewhere in the code
- Python driver handles: parameter passing, module loading, SLURM job script generation,
  result visualization
- C++ code compiled via CMake; `CMakeLists.txt` finds libMesh/PETSc via environment variables

## Output
- Velocity and pressure fields in ExodusII format (`.e`, ParaView-compatible)
- Console: per-Newton-iteration residual norm, total iterations, drag/lift coefficients C_D, C_L
