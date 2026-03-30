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
- **PETSc** (preferred for linear algebra, KSP/PC solvers, DM)
- **Trilinos** (alternative or complement: ML/MueLu for AMG, Belos for Krylov)
- FEM framework options (to be decided): **FEniCS/DOLFINx** (Python+C++, native PETSc) or **libMesh** (pure C++, supports PETSc+Trilinos)
- Do NOT implement linear solvers, preconditioners, or mesh routines from scratch
- Load libraries via `module load` on the cluster; use placeholders until module names are confirmed:
  ```
  module load petsc/X.Y trilinos/X.Y gcc/X.Y openmpi/X.Y
  ```

## Code Structure
- `src/` — C++ solver source files
- `scripts/` — Python driver (`run_simulation.py`) for setup, job submission, post-processing
- `meshes/` — geometry and mesh files (Gmsh `.geo` / `.msh`)
- `results/` — output (XDMF/HDF5 files for ParaView)

## Code Style
- All functions and classes must have concise docstrings/Doxygen comments
- Comments only where logic is non-obvious
- Python driver handles: parameter passing, module loading, job script generation, result visualization
- C++ code compiled via CMake; `CMakeLists.txt` must find PETSc/Trilinos via environment variables

## Output
- Velocity and pressure fields in XDMF/HDF5 format (ParaView-compatible)
- Console: iteration count, residual norms, drag/lift coefficients on cylinder
