# 2D Channel Flow Past a Cylinder

Time-dependent simulation of incompressible viscous flow in a channel past a
circular cylinder, using the finite element method (Taylor-Hood P2/P1 elements).
Starting from rest, the inlet velocity ramps up smoothly and recirculation
eddies develop behind the cylinder.

Based on the Schafer & Turek (1996) DFG benchmark geometry.

## Quick start

```bash
# 1. Build (requires libMesh + PETSc; see "Dependencies" below)
mkdir build && cd build
cmake -DLIBMESH_DIR=/path/to/libmesh ..
cmake --build . -j4

# 2. Generate mesh
python ../meshes/generate_mesh.py --output ../meshes/channel.msh

# 3. Run
mpirun -n 1 src/channel_flow --mesh ../meshes/channel.msh

# 4. View results in ParaView
#    Open results/channel_flow.e
```

The solver writes ExodusII snapshots to `results/channel_flow.e` as a time
series.  Open this file in ParaView to animate the eddy development.

## Problem description

```
  u = parabolic                     stress-free outlet
      inlet         +-----------+
    ───────────>    |     o     | ───────────>
                    |  cylinder |
      (BID 1)       +-----------+     (BID 2)
                      no-slip walls    (BID 3)
                      no-slip cylinder (BID 4)

  Channel: [0, 2.2] x [0, 0.41] m
  Cylinder: centre (0.2, 0.2), radius 0.05 m
```

The inlet velocity follows a parabolic profile that ramps up smoothly from
zero over `T_RAMP` seconds using a sin^2 function, allowing visualisation of
the transient eddy development.

## Parameters

All numerical parameters are centralised in `src/params.h`.  The most
important ones to tweak:

### Reynolds number

The Reynolds number is controlled by `U_MAX` (peak inlet velocity):

```
Re = U_MEAN * D / NU = (2/3 * U_MAX) * 0.1 / 0.001
```

| U_MAX (m/s) | Re  | Flow regime |
|-------------|-----|-------------|
| 0.075       | 5   | Symmetric eddies, very short |
| 0.15        | 10  | Symmetric eddies, moderate length (default) |
| 0.3         | 20  | Longer eddies, still steady |
| 1.5         | 100 | Vortex shedding (von Karman street) |

To change Re, edit `U_MAX` in `src/params.h` and rebuild.

### Cylinder position

The cylinder centre is set in two places:

1. **`src/params.h`**: `CYL_X`, `CYL_Y`, `CYL_RADIUS` (used by the solver
   for boundary conditions and drag/lift computation).

2. **`meshes/generate_mesh.py`**: `_XC`, `_YC`, `_R` (used for mesh
   generation).

Both must be updated consistently.  After changing, regenerate the mesh and
rebuild the solver.

### Time stepping

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DT` | 0.025 s | Time step |
| `T_FINAL` | 8.0 s | Final simulation time |
| `T_RAMP` | 2.0 s | Inlet velocity ramp-up duration |
| `OUTPUT_INTERVAL` | 10 | Write snapshot every N steps |
| `THETA` | 1.0 | 1.0 = backward Euler, 0.5 = Crank-Nicolson |

### Mesh resolution

Generate meshes with different resolutions using:

```bash
# Coarse (fast, for testing)
python meshes/generate_mesh.py --lc-far 0.1 --lc-cyl 0.03 --output meshes/channel_coarse.msh

# Medium (default)
python meshes/generate_mesh.py --lc-far 0.05 --lc-cyl 0.01 --output meshes/channel.msh

# Fine (for publication-quality results)
python meshes/generate_mesh.py --lc-far 0.02 --lc-cyl 0.003 --output meshes/channel_fine.msh
```

`--lc-far` controls the global element size; `--lc-cyl` controls the element
size near the cylinder surface.

## Visualisation

### ParaView (local)

1. Open `results/channel_flow.e` in ParaView.
2. Apply the "Warp By Vector" filter or colour by velocity magnitude.
3. Use the time controls (play button / time slider) to animate the eddy
   development.
4. For streamlines: Filters > Stream Tracer, seed along a vertical line
   upstream of the cylinder.

### ParaView from a remote cluster

If the simulation ran on a remote cluster (e.g. Oscar at Brown), there are
several ways to visualise the results:

**Option 1: Download the file**

```bash
# From your local machine:
scp user@ssh.ccv.brown.edu:/path/to/channel-flow-cylinder/results/channel_flow.e .
# Then open in local ParaView
```

The ExodusII `.e` file is typically 1-50 MB depending on mesh size and number
of time steps.

**Option 2: ParaView in client-server mode**

This avoids transferring the data.  On the cluster:

```bash
# Start ParaView server on the cluster
module load paraview
pvserver --server-port=11111
```

On your local machine, forward the port and connect:

```bash
ssh -L 11111:localhost:11111 user@ssh.ccv.brown.edu
```

Then in local ParaView: File > Connect > Add Server > `localhost:11111`.

**Option 3: VNC / Open OnDemand**

Many HPC centres (including Brown's Oscar) provide a web-based desktop via
Open OnDemand.  Launch a desktop session and run ParaView directly in the
remote desktop.  On Oscar: https://ood.ccv.brown.edu

### Useful ParaView tips

- Colour by `u` for x-velocity, `v` for y-velocity, or compute velocity
  magnitude with the Calculator filter: `sqrt(u^2 + v^2)`.
- To see the pressure field, colour by `p`.
- Use "Rescale to Data Range" at each time step for best contrast.
- Export animations: File > Save Animation (choose PNG sequence or AVI).

## Python driver

The driver script automates mesh generation and SLURM job submission:

```bash
# Dry run (show what would happen without executing)
python scripts/run_simulation.py --dry-run

# Generate mesh and submit to SLURM
python scripts/run_simulation.py --re 10 --dt 0.025 --t-final 8.0

# Use an existing mesh
python scripts/run_simulation.py --mesh-file meshes/channel.msh --ntasks 4
```

Run `python scripts/run_simulation.py --help` for all options.

## Dependencies

- **libMesh** (>= 1.7) with PETSc backend
- **PETSc** (>= 3.15) with Hypre (BoomerAMG)
- **Gmsh** Python package (`pip install gmsh`) for mesh generation
- **CMake** (>= 3.16)
- **MPI** runtime (OpenMPI or MPICH)

On Oscar (Brown CCV), these are available via spack:

```bash
source ~/data/spack/share/spack/setup-env.sh
spack load libmesh openmpi
```

## Testing

```bash
# C++ tests (from build directory)
cmake --build . --target generate_coarse_mesh   # one-time mesh generation
ctest --output-on-failure

# Python tests
python -m pytest tests/ -v
```

## Project structure

```
src/
  params.h                  Centralised numerical parameters
  channel_flow_system.h/cpp FEMSystem: variables, BCs, fieldsplit setup
  channel_flow_assembly.cpp Weak-form assembly (momentum, continuity, mass)
  main.cpp                  Time-stepping loop, ExodusII output
  drag_lift.h/cpp           Boundary stress integration for C_D, C_L
  petsc_utils.h             PETSc string formatting helper

meshes/
  generate_mesh.py          Gmsh mesh generation script

scripts/
  run_simulation.py         Python driver (mesh + SLURM + submission)

tests/
  test_spaces.cpp           FEM space initialisation
  test_stokes.cpp           Stokes linear solve
  test_output.cpp           ExodusII I/O + drag/lift
  test_timestep.cpp         Time-dependent stepping + ramp
  test_mesh.py              Mesh validity checks
  test_driver.py            Python driver tests
  test_validation.py        Schafer-Turek benchmark (Re=20 steady)
```

## Solver details

- **Spatial discretisation**: Taylor-Hood P2/P1 (inf-sup stable)
- **Time integration**: Backward Euler (theta=1.0), first-order
- **Nonlinear solver**: Newton with backtracking line search
- **Linear solver**: FGMRES with block Schur complement preconditioner
  - Velocity block: BoomerAMG (M/dt regularises the operator)
  - Pressure Schur complement: BoomerAMG on assembled Sp (SPD)
  - Lower-triangular Schur factorisation

## References

- Schafer, M. & Turek, S. (1996). "Benchmark computations of laminar flow
  around a cylinder." *Notes on Numerical Fluid Mechanics*, 52, 547-566.
- Elman, H., Silvester, D. & Wathen, A. (2014). *Finite Elements and Fast
  Iterative Solvers*, 2nd ed. Oxford University Press.
