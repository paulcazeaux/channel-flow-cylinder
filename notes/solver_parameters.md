# Solver Parameter Justification

## Finite Element Choice: Taylor-Hood P2/P1

Taylor-Hood elements (quadratic velocity, linear pressure) are the canonical choice
for incompressible Navier-Stokes:
- Satisfy the inf-sup (LBB) condition (Brezzi & Falk 1991)
- Second-order accurate in velocity, first-order in pressure
- No pressure stabilization (SUPG, PSPG) needed

**References:**
- Hood & Taylor (1974), "A numerical solution of the Navier-Stokes equations using the
  finite element technique", Computers & Fluids 1(1).
- Brezzi & Falk (1991), "Stability of higher-order Hood-Taylor methods", SIAM J. Numer. Anal.

## Mesh Sizing

- Global characteristic length h = 0.05 m
- Near cylinder h = 0.01 m (5× refinement)
- Rationale: boundary layer thickness δ ~ D/Re^0.5 ≈ 0.1/√5 ≈ 0.045 m at Re=5;
  h_cyl=0.01 gives ~4 elements across δ, sufficient for Re=5

**Reference:**
- Schafer & Turek (1996) recommend h ≈ 0.02 near cylinder for their benchmark meshes.

## Nonlinear Solver: Newton with Backtracking Line Search

- Newton's method gives quadratic convergence near the solution
- Line search (Armijo backtracking) ensures global convergence from a Stokes initial guess
- At Re=5 (strongly dominated by viscosity), Newton converges in 5–15 iterations

PETSc options:
```
-snes_type newtonls
-snes_linesearch_type bt      # backtracking (Armijo condition)
-snes_atol 1e-8
-snes_rtol 1e-10
-snes_max_it 50
```

**Reference:**
- Deuflhard (2004), "Newton Methods for Nonlinear Problems", Springer.
- Girault & Raviart (1986), "Finite Element Methods for Navier-Stokes Equations", Springer.

## Initialization: Stokes → Navier-Stokes

Solving the linearized Stokes problem first (advection term disabled) provides an
excellent initial guess for Newton's method:
- Stokes solution is unique and always exists
- Navier-Stokes at Re=5 is a small perturbation of Stokes; the two solutions are close
- Avoids Newton divergence that occurs from zero initial guess with nonzero advection

## Linear Solver: FGMRES + BoomerAMG

**FGMRES (Flexible GMRES, Saad 1993):**
- Required when preconditioner is variable (e.g., BoomerAMG uses inexact solves internally)
- Restart = 100: balances memory and convergence for ~10k–100k DOF systems
- Tolerances: rtol = 1e-10, max_it = 500

PETSc options:
```
-ksp_type fgmres
-ksp_gmres_restart 100
-ksp_rtol 1e-10
-ksp_max_it 500
```

**BoomerAMG (Hypre algebraic multigrid):**
- Optimal O(N) complexity for elliptic problems
- Effective for the viscosity-dominated velocity block at Re=5
- Standard settings sufficient; no tuning needed at low Re

PETSc options:
```
-pc_type hypre
-pc_hypre_type boomeramg
```

**References:**
- Saad (1993), "A flexible inner-outer preconditioned GMRES algorithm", SIAM J. Sci. Comput.
- Henson & Yang (2002), "BoomerAMG: a parallel algebraic multigrid solver", Appl. Num. Math.
- Benzi, Golub & Liesen (2005), "Numerical solution of saddle point problems", Acta Numerica.

## Drag/Lift Computation

Drag and lift coefficients are computed by boundary integration of the stress tensor
over the cylinder surface:

  F_D = ∫_Γ (ν ∂u/∂n - p·n_x) dΓ
  F_L = ∫_Γ (ν ∂v/∂n - p·n_y) dΓ

Normalization (Schafer-Turek):
  C_D = 2 F_D / (U_mean² D)
  C_L = 2 F_L / (U_mean² D)

where D = 2r = 0.1 m (cylinder diameter), U_mean = mean inlet velocity.
