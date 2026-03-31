# Implementation Plan — Time-Dependent Channel Flow Past a Cylinder

## Goal
Visualise recirculation eddies developing from rest in 2D channel flow past a
cylinder at Re = 5–10.  Produce an ExodusII time series for ParaView animation.

## Starting Point
The `main` branch has a complete, validated **steady-state** solver:
- Taylor-Hood P2/P1 FEMSystem with Dirichlet BCs
- Newton + FGMRES + fieldsplit Schur preconditioner (ILU velocity, AMG pressure)
- ExodusII output and drag/lift post-processing
- Validated at Re = 20 against Schafer-Turek DFG 2D-1: C_D = 5.571

## Key Changes from Steady to Time-Dependent

| Aspect | Steady (main branch) | Time-dependent (this branch) |
|--------|---------------------|------------------------------|
| Time solver | `SteadySolver` | `EulerSolver` or `Euler2Solver` (BDF2) |
| Mass matrix | Not assembled | `mass_residual()` override needed |
| Weak form | ν∫∇u·∇w + (u·∇)u·w − p∇·w = 0 | M ∂u/∂t + ν∫∇u·∇w + (u·∇)u·w − p∇·w = 0 |
| Velocity preconditioner | ILU(1) (Oseen is non-symmetric) | BoomerAMG (M/dt regularises → diag-dominant) |
| Initial condition | Stokes solve | Zero velocity (start from rest) |
| Output | Single ExodusII file | Time series: snapshots every N steps |
| Drag/lift | Single value | Time history C_D(t), C_L(t) printed per step |

## Phases

| # | Phase | Status |
|---|-------|--------|
| 1 | Mass matrix assembly (`mass_residual`) | ✅ Done |
| 2 | Time-stepping in `main.cpp` (backward Euler loop, dt, T_final) | ✅ Done |
| 3 | Preconditioner update (AMG for velocity block) | ✅ Done |
| 4 | Time-series ExodusII output (snapshots every N steps) | ✅ Done |
| 5 | Per-step drag/lift output C_D(t), C_L(t) | ✅ Done |
| 6 | Testing and validation | ✅ Done |
| 7 | Python driver update (time-dependent parameters) | ⬜ Todo |

## Phase Details

### Phase 1 — Mass matrix assembly

Add `element_mass_residual()` override to `ChannelFlowSystem`.  libMesh's
`UnsteadySolver` calls this to evaluate M · ∂u/∂t.

The mass residual for incompressible NS is:
```
∫ ρ (∂u/∂t) · w  dx     (momentum equations only; no time derivative in continuity)
```

Implementation: override `bool mass_residual(bool request_jacobian, DiffContext&)`.
Only velocity variables (u, v) contribute; pressure has no mass term.

**Files**: `channel_flow_system.h`, `channel_flow_assembly.cpp`

### Phase 2 — Time-stepping loop

Replace `SteadySolver` with libMesh's `Euler2Solver` (Crank-Nicolson/BDF2) or
`EulerSolver` (backward Euler) in `main.cpp`.

Key parameters (add to `params.h`):
```
DT        = 0.025      // time step [s]
T_FINAL   = 8.0        // final time [s] (several convective times L/U ≈ 22s at Re=5)
THETA     = 1.0        // implicit Euler (θ=1) or Crank-Nicolson (θ=0.5)
```

At Re = 5–10 the flow reaches steady state after ~3–5 convective times.
Convective time scale: L/U_mean ≈ 2.2/0.1 = 22 s at Re = 5.  T_final = 8 s
is sufficient to see eddies fully developed.

The time loop: `for (t = 0; t < T_FINAL; t += DT) { sys.solve(); advance(); }`.
libMesh's `UnsteadySolver` manages old solution vectors automatically.

Update `params.h` with U_MAX for Re = 5–10.

**Files**: `main.cpp`, `params.h`

### Phase 3 — Preconditioner update

With the mass matrix M/dt in the velocity block, the operator becomes:
```
A = M/dt + νK + N(u)
```
At small dt and moderate Re, M/dt dominates → the block is diagonally dominant
and nearly SPD.  BoomerAMG now works well for the velocity sub-block.

Change velocity sub-PC options:
```
-fieldsplit_velocity_pc_type  hypre
-fieldsplit_velocity_pc_hypre_type  boomeramg
```

Keep pressure sub-PC as BoomerAMG on assembled Sp (unchanged).

**Files**: `main.cpp`, `tests/test_stokes.cpp` (if updated)

### Phase 4 — Time-series ExodusII output

Write snapshots at regular intervals using `ExodusII_IO::write_timestep()`:
```cpp
ExodusII_IO exo(mesh);
exo.write_timestep(filename, es, timestep_number, time);
```

All time steps are appended to a single `.e` file, which ParaView can load
as an animation.

Add parameter: `OUTPUT_INTERVAL = 10` (write every 10 time steps).

**Files**: `main.cpp`, `params.h`

### Phase 5 — Per-step drag/lift

Call `compute_drag_lift()` at each output step and print C_D(t), C_L(t).
Existing `drag_lift.cpp` works unchanged — it reads from `current_local_solution`.

Optionally write a `drag_lift.csv` file for plotting.

**Files**: `main.cpp`

### Phase 6 — Testing and validation

- **test_timestep.cpp**: Run 5 time steps on coarse mesh, verify:
  - Newton converges at each step (≤ 10 iterations)
  - Solution changes between steps (transient is non-trivial)
  - C_D(t) is finite at final step
  - ExodusII output file has multiple time steps

- Update existing tests as needed (test_stokes may need adjustment for
  mass_residual changes).

**Files**: `tests/test_timestep.cpp`, `tests/CMakeLists.txt`

### Phase 7 — Python driver update

Update `scripts/run_simulation.py` to accept time-dependent parameters:
`--dt`, `--t-final`, `--output-interval`, `--re` (adjusts U_MAX).

**Files**: `scripts/run_simulation.py`, `tests/test_driver.py`

## Parameters (params.h additions)

```cpp
// Time-dependent solver
constexpr double DT           = 0.025;  // time step [s]
constexpr double T_FINAL      = 8.0;    // final simulation time [s]
constexpr int    OUTPUT_INTERVAL = 10;  // write ExodusII every N steps
constexpr double THETA        = 1.0;    // 1.0 = backward Euler, 0.5 = Crank-Nicolson

// Re = 5–10: adjust U_MAX
// Re = U_MEAN * D / NU = (2/3 U_MAX) * 0.1 / 0.001
// Re = 5:  U_MAX = 0.075 m/s
// Re = 10: U_MAX = 0.15  m/s
```

## Physical Expectations

At Re = 5–10:
- Flow starts from rest (u = v = 0 everywhere, with inlet BC applied)
- Boundary layer develops on cylinder, flow separates
- Symmetric recirculation eddies form in the wake
- Eddies grow and reach steady-state length
- Eddy length ≈ 0.1–0.3 D at Re = 5, ≈ 0.5–1.0 D at Re = 10
- C_D decreases from high initial value to steady state
- C_L ≈ 0 (symmetric flow, no vortex shedding below Re ≈ 47)

## File Map (changes from main branch highlighted)

```
src/
├── params.h                     # + DT, T_FINAL, OUTPUT_INTERVAL, THETA, Re adjust
├── channel_flow_system.h        # + mass_residual declaration
├── channel_flow_assembly.cpp    # + mass_residual implementation
├── main.cpp                     # rewritten: time loop replaces steady solve
├── drag_lift.h / drag_lift.cpp  # unchanged
├── petsc_utils.h                # unchanged
└── CMakeLists.txt               # unchanged

tests/
├── test_timestep.cpp            # new: Phase 6 time-stepping test
├── test_spaces.cpp              # unchanged
├── test_stokes.cpp              # may need minor update
├── test_output.cpp              # may need minor update
└── CMakeLists.txt               # + test_timestep

scripts/
└── run_simulation.py            # + time-dependent CLI args
```
