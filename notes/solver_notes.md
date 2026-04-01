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

---

## DG Solver: Pressure Schur Complement with Cahouet-Chabard

### Motivation: why not monolithic fieldsplit for DG

The CG solver uses PETSc fieldsplit with a Schur complement *preconditioner*:
FGMRES on the full (u,p) system, preconditioned by an approximate block
factorization.  This works well because PETSc handles the block scattering.

For DG, libMesh's FEMSystem does not scatter the DG neighbor residual and
4-block Jacobian (elem-elem, elem-neighbor, neighbor-elem, neighbor-neighbor)
into the global system.  This was verified experimentally in Session 17:
`dg_terms_are_active()` is never true, and even after manual
`neighbor_side_fe_reinit()`, the neighbor contributions are lost.

Switching to deal.II resolves the assembly issue (MeshWorker::mesh_loop
handles DG neighbor assembly natively).  We also switch the solver strategy:
instead of a monolithic fieldsplit preconditioner, we explicitly solve the
pressure Schur complement.

### Block structure

After DG assembly with SIPG viscous + Lax-Friedrichs advection + pressure
coupling + pressure-jump stabilization, the Newton system has block form:

```
[ A   B ] [ δu ]   [ f ]
[ B^T C ] [ δp ] = [ g ]
```

where:
- A = M_V/dt + θ(νK + N(u_k) + N'(u_k))  — velocity block
  - M_V: DG velocity mass matrix (block-diagonal for DG)
  - K: DG viscous stiffness (SIPG: volume + face penalty + symmetry)
  - N: advection operator (Lax-Friedrichs flux)
  - N': advection Jacobian
  - θ = 0.5 for Crank-Nicolson
- B: pressure gradient operator (volume ∫p∇·w + face {p}[[w·n]])
- B^T: velocity divergence operator (volume -∫q∇·u + face -{q}[[u·n]])
- C: pressure-jump stabilization γ_p·h·∫_F [[p]][[q]] dS

### Schur complement elimination

Eliminate velocity from the second block row:

```
S·δp = g - B^T·A^{-1}·f       where S = C - B^T·A^{-1}·B
```

Then recover velocity:

```
A·δu = f - B·δp
```

S is the pressure Schur complement.  It is never formed explicitly — each
application S·v requires one A-solve.

### Spectral properties of S

For Stokes (N=0) with backward Euler (θ=1, A = M_V/dt + νK):

- Small dt regime (M_V/dt dominates A):
    A ≈ M_V/dt  →  S ≈ C - dt·B^T·M_V^{-1}·B = C - dt·L_p
  S is dominated by the discrete pressure Laplacian L_p = B^T·M_V^{-1}·B.

- Large dt / steady regime (νK dominates A):
    A ≈ νK  →  S ≈ C - (1/ν)·B^T·K^{-1}·B ≈ C - (1/ν)·M_p
  S is dominated by the pressure mass matrix M_p (from continuous inf-sup).

The full Schur complement interpolates between these regimes.  For NS,
advection modifies the picture but does not change the dominant spectral
behavior (the Oseen operator M/dt + νK + N has similar block structure).

### Cahouet-Chabard preconditioner (Guermond et al.)

The preconditioner for S^{-1} that captures both spectral limits:

```
P_CC = ν · M_p^{-1} + (1/dt) · L_p^{-1}
```

References:
- Cahouet & Chabard (1988), "Some fast 3D finite element solvers for the
  generalized Stokes problem"
- Guermond, Minev & Shen (2006), "An overview of projection methods for
  incompressible flows" (§4.2)

### Why DG makes Cahouet-Chabard particularly efficient

For CG (Taylor-Hood P2/P1), both M_V and M_p are globally coupled sparse
matrices.  Computing M_V^{-1} in L_p = B^T·M_V^{-1}·B is expensive —
typically approximated by mass lumping or a few CG iterations.

For DG, mass matrices are **block-diagonal** (one dense block per element,
no inter-element coupling):

```
M_V = diag(M_V^{e1}, M_V^{e2}, ..., M_V^{eN})
M_p = diag(M_p^{e1}, M_p^{e2}, ..., M_p^{eN})
```

Block sizes for P2 on triangles: 6×6 (velocity per component), 6×6 (pressure).

Consequences:
1. M_V^{-1} is **exact**: invert each small dense block once at setup.
   Cost: O(n_elem · k^6) where k is polynomial order.  Negligible.
2. M_p^{-1} is **exact**: same.  Each P_CC application of the M_p^{-1}
   term is a single element-local dense matrix-vector product.
3. L_p = B^T·M_V^{-1}·B can be **assembled explicitly** as a sparse
   matrix at setup time (or reassembled when A changes for NS).
   This is practical because M_V^{-1} is cheap.
4. L_p^{-1} is approximated by 2 BoomerAMG V-cycles.  L_p is SPD and
   Laplacian-like → AMG is near-optimal.

### Pressure-jump stabilization interaction

The stabilization C = γ_p·h·[[p]][[q]] contributes to S but is small
(γ_p = 0.05).  Two strategies:

Option A: Ignore C in the preconditioner.  P_CC as defined above.
  Rationale: C is O(γ_p·h) while L_p is O(1/h) — C is asymptotically
  negligible on fine meshes.

Option B: Fold C into L_p:  P_CC = ν·M_p^{-1} + (1/dt)·(L_p + C)^{-1}
  Rationale: C has the same sparsity as L_p (face-coupled pressure DOFs),
  so L_p + C is still SPD and AMG-friendly.  Marginally better conditioning.

Start with Option A; switch to B if outer iteration counts exceed ~30.

### Full solve algorithm (one time step)

```
Given u_n, solve for (u_{n+1}, p_{n+1}):

  Pre-compute (once per time step or once overall for Stokes):
    - Invert M_V element-locally: store M_V^{-1} blocks
    - Invert M_p element-locally: store M_p^{-1} blocks
    - Assemble L_p = B^T · M_V^{-1} · B as sparse matrix
    - Set up AMG on L_p (BoomerAMG, strong threshold 0.7)

  For each Newton step (typically 2-3 for NS, 1 for Stokes):
    1. Assemble A, B, C, f, g at current (u_k, p_k)
    2. Set up velocity preconditioner on A (AMG or ILU)

    3. Velocity pre-solve:      A z = f              (AMG/ILU, iterative)
    4. Schur complement RHS:    r = g - B^T z         (matvec)
    5. Pressure solve:          S p = r               (GMRES + P_CC)
         Each GMRES iteration:
           - S·v:  w = B·v; solve A·y = w; return C·v - B^T·y
           - P_CC·v:  return ν·M_p^{-1}·v + (1/dt)·L_p^{-1}·v
                       (element-local + 2 AMG V-cycles)
    6. Velocity back-solve:     A u = f - B·p         (AMG/ILU, iterative)

    7. Update: u_{k+1} = u_k + δu,  p_{k+1} = p_k + δp
```

### Expected iteration counts

| Component | Method | Expected iterations |
|-----------|--------|-------------------|
| Newton (outer) | Newton-Raphson | 1 (Stokes), 2-3 (NS Re=100) |
| Pressure Schur (step 5) | GMRES + P_CC | 10-30 |
| Each A-solve inside step 5 | AMG or ILU | 20-50 (inexact is fine) |
| Velocity pre/back-solve (steps 3, 6) | AMG or ILU | 20-50 |

The pressure solve (step 5) dominates: each of its ~20 GMRES iterations
requires one A-solve (~30 inner iterations).  Total linear work per Newton
step: ~20 × 30 = 600 effective Krylov iterations.  Comparable to the
monolithic approach (50 outer FGMRES × ~10 inner work = 500).

### Comparison with CG monolithic fieldsplit

| Aspect | CG monolithic | DG Schur complement |
|--------|--------------|-------------------|
| Outer solver | FGMRES on full (u,p) | GMRES on pressure S |
| Velocity sub-PC | ILU(1) (1 sweep) | AMG or ILU (iterative) |
| Pressure sub-PC | AMG on selfp Sp | Cahouet-Chabard (AMG on L_p) |
| Mass matrix inverse | Implicit (in ILU) | Exact (DG block-diagonal) |
| Schur complement | Approximate (selfp) | Exact (nested Krylov) |
| Inter-element coupling | Shared DOFs | Face fluxes only |

The DG approach has more DOFs (~5× for same mesh) but the block-diagonal
mass matrices and exact Schur complement may compensate through lower
iteration counts and better parallel scalability.

## IMEX Time Integration: Explicit Advection, Implicit Stokes

### Motivation

With fully implicit time stepping, the velocity block
A = M_V/dt + θ(νK + N(u)) is non-symmetric due to the advection operator
N(u).  At Re=100 the cell Peclet number Pe_h = |u|h/ν reaches 50–140,
and AMG fails on this non-symmetric block (Session 16).  Workarounds
(ILU, block-Jacobi+GMRES) are robust but not scalable.

IMEX treats advection explicitly and viscous+pressure implicitly.  The
implicit velocity block becomes:

    A = M_V/dt + θνK_DG       (symmetric positive definite)

AMG is near-optimal for this operator at any Reynolds number.  Furthermore,
A is constant across a time step (no dependence on the current velocity
iterate), so the AMG hierarchy is set up once and reused.

### IMEX scheme: high-order for larger CFL

We use a high-order IMEX scheme to maximize the stable time step.  Higher
order in the explicit part gives a larger CFL constant for the same formal
accuracy.

**IMEX-RK (Ascher, Ruuth, Spiteri 1997; Kennedy & Carpenter 2003):**

IMEX Runge-Kutta methods pair an explicit RK tableau (for advection) with
a diagonally-implicit RK (DIRK) tableau (for diffusion + pressure).  At
each stage s, the implicit solve is:

    (M_V + a_ss · dt · νK_DG) u^{(s)} + a_ss · dt · B p^{(s)}
        = M_V u_n - dt Σ_{j<s} ã_sj N(u^{(j)})      ← explicit advection
                   - dt Σ_{j<s} a_sj (νK u^{(j)} + B p^{(j)})

    B^T u^{(s)} + C p^{(s)} = 0

Each stage solves the *same* operator A = M_V + a_ss·dt·νK_DG (for SDIRK
methods where a_ss is constant across stages).

Candidate schemes:

| Scheme | Order | Stages | CFL (advection) | a_ss |
|--------|-------|--------|-----------------|------|
| IMEX-SSP2(2,2,2) | 2 | 2 | ~1.0 | 0.5 |
| IMEX-ARK3(4,4,3) | 3 | 4 | ~1.7 | 0.435867 |
| IMEX-ARK4(6,6,4) | 4 | 6 | ~2.0 | — |

Higher order → more stages per step, but each stage uses the same A
operator with the same AMG setup.  The larger CFL constant partially
offsets the smaller dt requirement.

For our problem: with 3rd-order IMEX and CFL ≈ 1.7:
    dt ≤ 1.7 · h_min / |u_max| ≈ 1.7 · 0.01 / 1.5 ≈ 0.011 s

Compared to CFL ≈ 1.0 for first-order:
    dt ≤ 1.0 · 0.01 / 1.5 ≈ 0.007 s

The 3rd-order scheme allows ~60% larger dt.

### IMEX-ARK3(4,4,3) — Kennedy & Carpenter (2003)

This is a popular choice: 3rd order, 4 stages, L-stable implicit part,
SSP explicit part.  The DIRK coefficient a_ss = (1+√2/2)^{-1} ≈ 0.4359
is the same for all diagonal entries (SDIRK), so one AMG setup suffices.

Explicit Butcher tableau (advection):

    0     |
    c2    | ã21
    c3    | ã31  ã32
    1     | ã41  ã42  ã43
    ------+-----------------
          | b1   b2   b3   b4

Implicit Butcher tableau (diffusion + pressure):

    0     | 0
    c2    | a21   a22
    c3    | a31   a32   a33
    1     | a41   a42   a43   a44
    ------+-------------------------
          | b1    b2    b3    b4

With a22 = a33 = a44 = γ (SDIRK property).
The implicit operator at every stage is:

    A = M_V + γ·dt·νK_DG       (constant across all stages and time steps)

### Impact on Schur complement solver

With IMEX, each implicit stage solves a linear Stokes-like system (no
Newton iteration needed):

    [ A   γ·dt·B ] [ u^{(s)} ]   [ RHS_u ]
    [ B^T  C      ] [ p^{(s)} ] = [ RHS_p ]

where A = M_V + γ·dt·νK_DG is SPD.  The Schur complement is:

    S = C - B^T A^{-1} B

and the Cahouet-Chabard preconditioner is:

    P_CC = ν · M_p^{-1} + 1/(γ·dt) · L_p^{-1}

Everything computed at setup time:
- M_V^{-1}, M_p^{-1}: element-local dense inverses (DG block-diagonal)
- L_p = B^T M_V^{-1} B: assembled once as sparse matrix
- AMG on L_p: set up once
- AMG on A: set up once (A is constant)

Per-stage cost:
1. Evaluate explicit advection N(u^{(j)}) at previous stages (matvec)
2. Form RHS from Butcher tableau coefficients
3. One Schur complement solve (GMRES + P_CC, ~15–25 iterations)
4. Each Schur iteration: one A-solve (AMG, ~5–10 iterations)
5. Velocity back-solve (one A-solve)

Per time step: 4 stages × (1 Schur solve + 2 A-solves) = ~4 × 30 = 120 A-solves.

### Cost comparison: fully implicit vs IMEX-ARK3

Fully implicit (dt=0.02, 400 steps, 2.5 Newton/step):
- Per step: 2.5 × (20 Schur iters × 30 inner + 2 × 30) ≈ 1650 A-equivalent
- Total: 400 × 1650 = 660,000 A-equivalent
- A-solve: ~30 iters (ILU, non-symmetric) each

IMEX-ARK3 (dt=0.011, 727 steps, 4 stages/step):
- Per step: 4 × (20 Schur iters × 8 inner + 2 × 8) ≈ 704 A-equivalent
- Total: 727 × 704 = 512,000 A-equivalent
- A-solve: ~8 iters (AMG, SPD) each — much cheaper per iteration

IMEX is ~22% fewer total A-equivalents, and each A-equivalent is
cheaper (AMG vs ILU, SPD vs non-symmetric).  Net speedup: ~2-3×.
Additionally, AMG parallelizes well while ILU does not.

### Stability considerations

- The explicit advection is subject to CFL.  For vortex shedding at
  Re=100, the maximum velocity |u_max| ≈ 1.5 m/s occurs at the inlet
  and near the cylinder.  With h_min ≈ 0.01 near the cylinder:
  dt_CFL ≈ 1.7 × 0.01 / 1.5 ≈ 0.011 s.

- This is comparable to our current dt=0.02 (only ~2× smaller).
  The temporal accuracy at 3rd order with dt=0.011 is better than
  2nd order (Crank-Nicolson) with dt=0.02.

- The implicit SDIRK part is L-stable: no stability constraint from
  viscous or pressure terms.

- For adaptive time stepping: monitor CFL and adjust dt.  The AMG
  hierarchy on A must be rebuilt if γ·dt·ν changes significantly.
