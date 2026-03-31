# Solver and Preconditioner Notes

Summary of experiments, reasoning, and choices made during development.

## Time Integrator

### Backward Euler (θ=1, EulerSolver)

- First-order accurate, L-stable, unconditionally stable
- Requires small dt for temporal accuracy (dt ≈ 0.005 at Re=100)
- Newton converges quadratically: 2 iterations per step
- Heavy numerical damping can suppress physical oscillations

### Crank-Nicolson (θ=0.5, EulerSolver)

- Second-order accurate, A-stable (not L-stable)
- Allows 4× larger dt for equivalent accuracy (dt ≈ 0.02 at Re=100)
- Newton converges quadratically when the Jacobian is assembled correctly
  (see "Jacobian correctness" below): 2 iterations per step

**Chosen: Crank-Nicolson (θ=0.5)** for second-order accuracy.

### Euler2Solver (θ=0.5)

Evaluated and rejected. Euler2Solver evaluates θ·F(u_new) + (1-θ)·F(u_old)
rather than F(θ·u_new + (1-θ)·u_old). The Jacobian only captures the u_new
contribution, giving inherently **linear Newton convergence** with factor
(1-θ) = 0.5 — 10 iterations per step. Not fixable by better Jacobian assembly.

### BDF2

Not available in libMesh 1.7. Would require a custom implementation storing
two previous solution vectors. The EulerSolver with θ=0.5 provides second-order
accuracy without this complication.

## Jacobian Correctness for Crank-Nicolson

Two bugs were found in the original assembly that degraded Newton convergence
from quadratic (2 iterations) to linear (10 iterations).

### Bug 1: Missing θ factor in element_time_derivative

libMesh's EulerSolver sets `context.elem_solution_derivative = θ` before
calling `element_time_derivative`. The DiffContext documentation says:
"Corresponding Jacobian contributions should be multiplied by this amount."

This is the chain-rule factor ∂u_θ/∂u = θ, needed because the residual
F(u_θ) depends on u through the blended state u_θ = θu + (1-θ)u_old.

**Fix**: multiply all Jacobian entries in `element_time_derivative` by
`context.get_elem_solution_derivative()`.

Note: `element_constraint` does NOT need this scaling — EulerSolver calls
it at u_new (not u_θ) with `elem_solution_derivative = 1` (see
euler_solver.C line 179).

### Bug 2: Wrong function and missing 1/dt in mass_residual

The mass residual must read the **rate** u̇ = (u_new − u_old)/dt, not the
blended state u_θ:

- **Wrong**: `c.interior_value(var, qp)` → reads u_θ from elem_solution
- **Correct**: `c.interior_rate(var, qp)` → reads u̇ from elem_solution_rate

The mass Jacobian must include the factor ∂u̇/∂u = 1/dt:

- **Wrong**: M (bare mass matrix)
- **Correct**: `c.get_elem_solution_rate_derivative() * M` = M/dt

Found by reading the libMesh source (euler_solver.C) and the
`fem_system_ex1/naviersystem.C` example, which uses `c.interior_rate()`.

## Block Preconditioner: Fieldsplit Schur Complement

The incompressible Navier-Stokes saddle-point system is preconditioned using
PETSc's PCFIELDSPLIT with a Schur complement factorisation:

```
P⁻¹ ≈ [ A   0 ]⁻¹
       [ B   S̃ ]
```

where S̃ = B diag(A)⁻¹ Bᵀ (the "selfp" Schur approximation, SPD).

**Lower-triangular** factorisation (not diagonal): captures the B·A⁻¹
coupling term, reducing FGMRES iterations from 176 to 35 (tested in
Session 9 at Re=20, steady).

**SELFP Schur approximation**: must be set programmatically via
`PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL)`.
The PETSc option string `-pc_fieldsplit_schur_pre_type selfp` is
silently ignored (parsing bug or limitation).

## Velocity Sub-PC: ILU vs BoomerAMG

### When AMG works (low Re, time-dependent)

At low Re (≤ 10) with time stepping, M/dt dominates the velocity block:

```
A = M/dt + νK + N(u)
M/dt ≈ 50,  νK ≈ 0.1,  N(u) ≈ 1   (at Re=10, dt=0.025)
```

The block is diagonally dominant and nearly SPD → BoomerAMG works well.

**AMG tuning** (following Guermond): strong threshold 0.1 for the
mass-dominated velocity block, 0.7 for the Laplacian-like pressure
Schur complement. This alone reduced FGMRES by 32%.

Chebyshev smoother was tested but degraded convergence on the
non-symmetric velocity block. Default hybrid symmetric Gauss-Seidel
is better.

Fixed V-cycle count (2 or 4 V-cycles via Richardson KSP) was too weak
as a sub-PC inside the Schur fieldsplit — FGMRES hit the 200-iteration
limit. The sub-solves need higher accuracy than in Guermond's context
(where AMG is used as a standalone solver, not nested inside a
block preconditioner).

### When AMG fails (high Re)

At Re=100 with dt=0.02, the cell Péclet number is:

```
Pe_h = U·h/ν = 1.0·0.05/0.001 = 50  (fine mesh)
Pe_h = 1.0·0.1/0.001 = 100           (coarse mesh)
```

The advection Jacobian N(u) dominates: the velocity block is strongly
non-symmetric. BoomerAMG (Falgout coarsening, symmetric GS smoother) is
designed for SPD systems and stalls: FGMRES hits 200 iterations, Newton
fails.

### When ILU works (high Re, steady, or time-dependent)

ILU(k) is algebraically robust for non-symmetric systems. It doesn't
require any spectral properties of the matrix.

Tested ILU(1) and ILU(2):

| Sub-PC | FGMRES/step (Re=100, fine mesh) | Scalability |
|--------|--------------------------------|-------------|
| ILU(1) | ~50 | Poor (fill grows with mesh) |
| ILU(2) | ~40 | Worse (more fill) |
| AMG | Diverged | Optimal when applicable |

ILU is the safe choice for production runs at Re ≥ 20. AMG is preferable
at low Re for its optimal (mesh-independent) complexity.

**Chosen**: ILU(1) for Re=100, BoomerAMG (threshold 0.1) for Re ≤ 10.

## Pressure Sub-PC

BoomerAMG with strong threshold 0.7 throughout. The pressure Schur
complement S̃ = B diag(A)⁻¹ Bᵀ is SPD and Laplacian-like — the ideal
AMG target. Never had convergence issues with this choice.

## FGMRES (Outer Krylov Solver)

- **FGMRES** (not GMRES): required because the preconditioner (AMG or
  ILU inside fieldsplit) is not exactly constant between iterations.
- Restart = 100, max iterations = 200.
- Linear tolerance: set adaptively by libMesh's inexact Newton
  (Eisenstat-Walker), not overridden. Overriding with a tight tolerance
  (e.g. -ksp_rtol 1e-10) forces unnecessary iterations.

## Summary of Choices by Reynolds Number

| Re | dt | Time integrator | Velocity PC | Pressure PC | Newton/step |
|----|------|----------------|-------------|-------------|-------------|
| 5-10 | 0.025 | CN (θ=0.5) | AMG (θ=0.1) | AMG (θ=0.7) | 2 |
| 20 | 0.02 | CN (θ=0.5) | ILU(1) | AMG (θ=0.7) | 2 |
| 100 | 0.02 | CN (θ=0.5) | ILU(1) | AMG (θ=0.7) | 2–3 |
| 100 (steady) | — | Newton | ILU(1) | AMG (θ=0.7) | 4 |
