# Implementation Plan — Vortex Shedding at Re = 100

## Goal
Visualise the von Karman vortex street behind a cylinder in channel flow at
Re = 100.  Produce an ExodusII time series for ParaView animation showing
periodic vortex shedding.

## Starting Point
The `time-dependent-ns` branch has a working time-dependent solver:
- Backward Euler time integration with sin² inlet ramp
- Mass matrix assembly (`mass_residual`)
- BoomerAMG for both velocity and pressure blocks
- ExodusII time-series output with drag/lift per step
- Validated at Re = 10 (symmetric eddies, steady state reached)

## Changes from Re = 10 to Re = 100

| Parameter | Re = 10 branch | Re = 100 (this branch) |
|-----------|---------------|----------------------|
| U_MAX | 0.15 m/s | 1.5 m/s |
| DT | 0.025 s | 0.005 s |
| T_RAMP | 2.0 s | 1.0 s |
| OUTPUT_INTERVAL | 10 | 20 |
| Expected flow | Symmetric eddies → steady | Vortex shedding → periodic |

## Physical Expectations (Re = 100)

- Vortex shedding begins after the initial transient (~1–2 s)
- Strouhal number St ≈ 0.2 → shedding period T_shed ≈ D/(St·U_MEAN) = 0.5 s
- C_D oscillates around a mean ≈ 3.2 with amplitude ≈ 0.01
- C_L oscillates with amplitude ≈ 1.0 and period ≈ 0.5 s
- Schafer-Turek 2D-2 reference values:
  - C_D_max ∈ [3.22, 3.24], C_L_max ∈ [0.99, 1.01]
  - St ∈ [0.295, 0.305] (based on cylinder diameter and U_MEAN)

## Phases

| # | Phase | Status |
|---|-------|--------|
| 1 | Parameter update (U_MAX, DT, T_RAMP, OUTPUT_INTERVAL) | ✅ Done |
| 2 | Run and verify vortex shedding develops | ⬜ Todo |
| 3 | Validate against Schafer-Turek 2D-2 benchmark | ⬜ Todo |
| 4 | Update README for vortex shedding case | ⬜ Todo |

## Solver Considerations

At Re = 100 with dt = 0.005:
- Cell Peclet: Pe_h = U·h/ν = 1.0·0.05/0.001 = 50 (convection-dominated)
- But M/dt ratio: M/(dt·νK) ≈ ρ/(dt·ν/h²) = 1/(0.005·0.001/0.0025) = 500

M/dt still dominates the velocity block → BoomerAMG should remain viable.
If the solver stalls, fall back to ILU for the velocity block.

The flow is now genuinely unsteady (not converging to steady state), so every
time step requires a real Newton solve.  Expect higher computational cost.

## Mesh Requirements

The fine mesh (lc_far=0.03, lc_cyl=0.005) should be sufficient.  At Re=100
the boundary layer is thinner (δ ~ D/√Re ≈ 0.01 m) so cylinder refinement
matters more.  If results are noisy, try lc_cyl=0.003.
