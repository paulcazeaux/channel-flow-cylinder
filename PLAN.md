# Implementation Plan — DG Navier-Stokes Solver (deal.II)

## Context

The existing CG Taylor-Hood (P2/P1) solver (libMesh) is fully validated against
Schafer-Turek benchmarks at Re=20 (steady) and Re=100 (vortex shedding).
We now implement a **Discontinuous Galerkin** formulation for the same problem
to enable high-order accuracy, local conservation, and future hp-adaptivity.

Starting point: the `vortex-shedding` branch (all 16 sessions of CG work).

### Why deal.II instead of libMesh

libMesh's `FEMSystem` framework does not support DG assembly out of the box.
The `DGFEMContext` class provides data structures (neighbor shape functions,
4-block Jacobians) but `FEMSystem::assembly()` does not scatter the neighbor
residual or DG Jacobian blocks into the global system. We verified this
experimentally: `dg_terms_are_active()` is never set by `FEMSystem`, and
even after manual `neighbor_side_fe_reinit()`, the neighbor contributions
are not assembled into the global matrix.

deal.II has first-class DG support:
- `FEFaceValues` / `FESubfaceValues` for face quadrature
- `MeshWorker::mesh_loop()` with `cell_worker`, `face_worker`, `boundary_worker`
  callbacks that automatically handle neighbor assembly
- `AffineConstraints` for hanging nodes and boundary conditions
- PETSc backend via `PETScWrappers` for the linear algebra
- Built-in support for `FE_DGQ` (tensor-product DG) and `FE_DGP` (simplex DG)

## Formulation: SIPG + Upwind DG, Equal-Order Q_k/Q_k

**SIPG** for viscous terms, **Lax-Friedrichs** for advection, **equal-order**
Q_k/Q_k with pressure-jump stabilization on interior faces.

Starting polynomial order: Q2/Q2. Can raise to Q3, Q4 by changing one parameter.

**Pressure null space:** Pin one pressure DOF to zero (simplest approach).

**Testing:** Same coarse mesh, 300s timeout.

## Time Integration: IMEX (Explicit Advection, Implicit Stokes)

Advection is treated explicitly; viscous + pressure implicitly.  This makes
the implicit velocity block A = M_V/(γ·dt) + νK_DG symmetric positive
definite, enabling optimal AMG solves regardless of Reynolds number.

### IMEX-ARK3 (Kennedy & Carpenter 2003)

3rd-order IMEX Runge-Kutta: 4 stages, L-stable SDIRK implicit part,
SSP explicit part.  CFL constant ≈ 1.7 (vs 1.0 for 1st-order IMEX).

For our problem: dt ≤ 1.7 · h_min / |u_max| ≈ 1.7 · 0.01 / 1.5 ≈ 0.011 s.

At each stage s, the implicit system is linear (no Newton needed):

```
[ A       γ·dt·B ] [ u^{(s)} ]   [ RHS_u(explicit advection) ]
[ (γ·dt)B^T  C   ] [ p^{(s)} ] = [ RHS_p                     ]
```

where A = M_V + γ·dt·νK_DG is SPD and constant across all stages and
time steps (SDIRK property: γ = a_ss is the same for all diagonal entries).
AMG is set up once and reused throughout.

Higher-order IMEX schemes (ARK4, ARK5) can further increase the CFL
constant at the cost of more stages per step.

## Linear Solver: Pressure Schur Complement + Cahouet-Chabard

### Schur complement solve (each IMEX stage)

```
1. Velocity pre-solve:    A z = RHS_u               (AMG, SPD)
2. Schur RHS:             r = RHS_p - (γ·dt)B^T z
3. Pressure solve:        S p = r                    (GMRES + P_CC)
4. Velocity back-solve:   A u = RHS_u - γ·dt·B p     (AMG, SPD)
```

S = C - (γ·dt)^2 B^T A^{-1} B is the Schur complement.

### Cahouet-Chabard preconditioner for S

```
P_CC = ν · M_p^{-1} + 1/(γ·dt) · L_p^{-1}
```

where L_p = B^T M_V^{-1} B (discrete pressure Laplacian).

DG advantages:
- M_p^{-1}, M_V^{-1}: exact, element-local (block-diagonal mass matrices)
- L_p: assembled explicitly once (cheap with exact M_V^{-1})
- L_p^{-1}: 2 AMG V-cycles (BoomerAMG, SPD Laplacian-like)
- A^{-1}: AMG (SPD, setup once, reused across all stages and time steps)

### Setup cost (once)

1. Invert M_V, M_p element-locally: store dense inverse blocks
2. Assemble L_p = B^T M_V^{-1} B as sparse matrix
3. Build AMG hierarchy on L_p
4. Build AMG hierarchy on A = M_V + γ·dt·νK_DG

### Per-stage cost

1. Evaluate explicit advection N(u^{(j)}) at previous stages (matvec)
2. Form RHS from Butcher tableau
3. Schur solve: ~15-25 GMRES iterations, each requiring one A-solve (~5-10 AMG iters)
4. Velocity back-solve: one A-solve

### Per-step cost (4 stages)

~4 × (20 Schur iters × 8 AMG iters + 2 × 8 AMG iters) ≈ 700 AMG-equivalents

## Phases

| # | Phase | Status |
|---|-------|--------|
| 0 | Install deal.II via spack | In progress |
| 1 | Mesh + DG space setup (FE_DGQ, DoFHandler) | |
| 2 | Block assembly: A, B, C, M_V, M_P via mesh_loop | |
| 3 | Explicit advection operator (Lax-Friedrichs DG flux) | |
| 4 | SIPG viscous + pressure face integrals (boundary + interior) | |
| 5 | Schur complement solver + Cahouet-Chabard preconditioner | |
| 6 | IMEX-ARK3 time stepper | |
| 7 | Drag/lift + VTK output | |
| 8 | Validation against Schafer-Turek | |

## Key deal.II Components

- `Triangulation<2>`: mesh (read from Gmsh via `GridIn`)
- `FE_DGQ<2>(k)`: tensor-product DG basis of order k
- `FESystem<2>`: blocked (u, v, p) or separate DoFHandlers per variable
- `DoFHandler<2>`: DG DOF management
- `MeshWorker::mesh_loop()`: automatic cell + face iteration with neighbor handling
- `FEValues<2>`, `FEFaceValues<2>`: quadrature point data
- `BlockSparseMatrix`, `BlockVector`: block linear algebra
- `SolverGMRES` + `LinearOperator` for Schur complement solve
- PETSc or deal.II native for AMG (TrilinosWrappers::PreconditionAMG or
  PETScWrappers if available)

## Files (new deal.II solver)

| File | Description |
|------|-------------|
| `src/dg_navier_stokes.h` | Solver class declaration |
| `src/dg_assembly.cc` | Assembly via mesh_loop (A, B, C, advection) |
| `src/dg_schur_solver.cc` | Schur complement + Cahouet-Chabard |
| `src/dg_time_stepper.cc` | IMEX-ARK3 time integration |
| `src/dg_params.h` | Centralized parameters |
| `src/main.cc` | Entry point |
| `src/CMakeLists.txt` | Build with deal.II |
| `tests/test_dg_spaces.cc` | DG DOF count verification |
| `tests/test_dg_stokes.cc` | Stokes solve + weak BC check |
