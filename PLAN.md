# Implementation Plan — Vortex Shedding at Re = 100

## Goal
Visualise the von Karman vortex street behind a cylinder in channel flow at
Re = 100.  Produce an ExodusII time series and MP4 movie showing periodic
vortex shedding.

## Phases

| # | Phase | Status |
|---|-------|--------|
| 1 | Parameter update (U_MAX=1.5, DT=0.02, θ=0.5) | ✅ Done |
| 2 | Crank-Nicolson time integrator + exact Newton Jacobian | ✅ Done |
| 3 | AMG tuning and preconditioner experiments | ✅ Done |
| 4 | Run and verify vortex shedding (refined mesh) | ✅ Done |
| 5 | Validate against Schafer-Turek 2D-2 benchmark | ✅ Done |
| 6 | Movie rendering (pvpython + ffmpeg) | ✅ Done |

## Validation (Refined Mesh, lc_far=0.02, lc_cyl=0.003)

| Quantity | Our result | Schafer-Turek 2D-2 |
|----------|-----------|-------------------|
| C_D_max | 3.223 | [3.22, 3.24] |
| C_L_max | 1.000 | [0.99, 1.01] |

## Solver Configuration

- Time integrator: EulerSolver θ=0.5 (Crank-Nicolson, 2nd order)
- dt = 0.02, T_final = 8.0, T_ramp = 1.0
- Newton: 2–3 iterations per step (quadratic convergence)
- Linear: FGMRES + fieldsplit Schur (lower, selfp)
  - Velocity: ILU(1) (AMG fails at Re=100, cell-Pe too high)
  - Pressure: BoomerAMG (strong threshold 0.7)

See `notes/solver_notes.md` for full details and experiment log.
