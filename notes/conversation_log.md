# Conversation Log

Chronological record of design decisions made across sessions.
Append new entries at the bottom; do not edit past entries.

---

## 2026-03-29 — Session 1: Project Setup

**Topics:** Project initialization, guidelines, implementation plan, scaffolding.

### Decisions made

**Architecture:**
- Solver: C++ with libMesh + PETSc (not DOLFINx, which is primarily Python-driven)
- Driver: Python for mesh generation, job submission, post-processing
- Rationale: user prefers pure C++ numerics with Python only as orchestration layer

**FEM choices (from literature):**
- Elements: Taylor-Hood P2/P1 — inf-sup stable, no stabilization needed at Re=5
- Weak form: standard Galerkin with skew-symmetric advection
- Nonlinear solver: PETSc SNES Newton + backtracking line search
- Linear solver: FGMRES (restart=100) + BoomerAMG (Hypre)
- Initialization: Stokes solution → Newton initial guess

**Domain:**
- Schafer-Turek 2D-1: channel [0,2.2]×[0,0.41], cylinder center (0.2,0.2) r=0.05
- Re=5 (primary), Re=20 for validation against published benchmark values

**Coding conventions established:**
- Source files ≤ 300 lines
- All numerical parameters in centralized `src/params.h`
- ExodusII (.e) output format
- PLAN.md as living phase tracker in repo
- `notes/` for literature justifications and benchmark targets

### Files created this session
- `CLAUDE.md` — project guidelines
- `PLAN.md` — implementation tracker
- `.gitignore`
- `notes/solver_parameters.md`
- `notes/benchmark_targets.md`
- `notes/conversation_log.md` (this file)

### Next steps (Phase 1)
- `meshes/channel.geo` — Gmsh parametric geometry
- `meshes/generate_mesh.py` — mesh generation script

---

## 2026-03-29 — Session 2: Mesh, Test Plan, CMake

**Topics:** Phase 1 (Gmsh mesh generation), test plan and testing guidelines,
Phase 5 (CMake build system).

### Decisions made

**Mesh generation:**
- `meshes/channel.geo` — parametric Gmsh geometry for the channel-cylinder domain
- `meshes/generate_mesh.py` — Python script using the Gmsh API; exposes `--lc-far`
  and `--lc-cyl` size parameters and `--output` path
- Physical boundary tags: 1 = inlet, 2 = outlet, 3 = walls (top+bottom), 4 = cylinder
- Production mesh: lc_far=0.05, lc_cyl=0.01; coarse test mesh: lc_far=0.1, lc_cyl=0.03

**Test framework:**
- C++ tests: standalone executables registered with CTest (`add_test`); mesh path
  passed via `TEST_MESH` environment variable so tests run without cluster access
- Python tests: `unittest` + `pytest`; coarse mesh generated as a CMake custom target
- Each phase has its own test file; 300-line limit applies to test files too
- Test plan added to `PLAN.md §Test Plan` with explicit pass/fail criteria per phase

**Build system (Phase 5 before Phase 2):**
- Top-level `CMakeLists.txt` uses `libmesh-config` to query all compile/link flags
- `src/CMakeLists.txt` defines `channel_flow_lib` (STATIC, linked by tests) and
  `channel_flow` executable; `SOLVER_SOURCES` list populated phase-by-phase
- `tests/CMakeLists.txt` defines `add_solver_test` macro and `generate_coarse_mesh`
  custom target; test executables uncommented as phases complete
- Phase 5 done before Phase 2 so the build infrastructure exists before solver code

### Files created this session
- `meshes/channel.geo`
- `meshes/generate_mesh.py`
- `tests/test_mesh.py` (Phase 1 test)
- `CMakeLists.txt` (top-level)
- `src/CMakeLists.txt`
- `tests/CMakeLists.txt`
- `src/params.h` — all numerical/physical constants
- `src/main.cpp` — Phase 5 build stub

### Next steps (Phase 2)
- `src/channel_flow_system.h` / `.cpp` — FEMSystem subclass with weak form
- `tests/test_spaces.cpp` — DOF space and Dirichlet BC verification

---

## 2026-03-29 — Session 3: Phase 2 FEMSystem + server bug fix

**Topics:** Phase 2 (Taylor-Hood FEMSystem assembly), merge of server-side fix.

### Decisions made

**FEMSystem design:**
- Subclasses `libMesh::FEMSystem`; overrides `init_data`, `init_context`,
  `element_time_derivative` (momentum), `element_constraint` (continuity)
- Variables added in `init_data`: u (SECOND/LAGRANGE), v (SECOND/LAGRANGE),
  p (FIRST/LAGRANGE); parent `FEMSystem::init_data()` called at end
- Dirichlet BCs applied via `DofMap::add_dirichlet_boundary`:
  - No-slip (walls + cylinder): `ZeroFunction` on both velocity components
  - Inlet u: `InletVelocityU` (anonymous-namespace `FunctionBase` subclass);
    `DirichletBoundary` clones the function internally so no lifetime issue
  - Inlet v: `ZeroFunction`
  - Outlet: natural Neumann, no action

**Weak form:**
- Standard convective form ∫(u·∇)u·w dx (equivalent to skew-symmetric for
  divergence-free u); sign convention: residual = 0 at solution
- Newton Jacobian assembled manually in both element methods (returned when
  `request_jacobian == true`); avoids finite-difference fallback
- Continuity residual sign: −∫q∇·u dx (consistent with off-diagonal Jacobian blocks)

**main.cpp:**
- Replaced Phase 5 stub with full solve loop: reads `--mesh` arg, constructs
  `SteadySolver` + `NewtonSolver` chain (tolerances from `params.h`), calls
  `sys.solve()`; Phase 3 will add PETSc KSP/PC flags, Phase 4 adds output

**Server-side fix merged:**
- Commit `1a48a00` (pushed from HPC cluster) fixed `test_no_degenerate_elements`
  in `tests/test_mesh.py`: `gmsh.model.mesh.getElements()` returns
  `(elementTypes, elementTags, nodeTags)` — the test was unpacking it as
  `(_, etypes, _)` (etypes held tag arrays, not type IDs), causing a TypeError.
  Fixed by unpacking correctly and passing tag arrays directly to
  `getElementQualities`.
- Pull was a clean fast-forward (no overlap with Phase 2 files).

### Files created/modified this session
- `src/channel_flow_system.h` — new
- `src/channel_flow_system.cpp` — new
- `tests/test_spaces.cpp` — new
- `src/main.cpp` — replaced stub
- `src/CMakeLists.txt` — enabled `channel_flow_system.cpp` in SOLVER_SOURCES
- `tests/CMakeLists.txt` — enabled `add_solver_test(test_spaces ...)`
- `PLAN.md` — Phase 2 marked done

### Next steps (Phase 3)
- Configure PETSc SNES/KSP: FGMRES (restart=100) + BoomerAMG via command-line flags
- Implement Stokes → Navier-Stokes initialization sequence in main.cpp
- `tests/test_stokes.cpp` — verify Stokes solve converges, residual < 1e-8

---

## 2026-03-30 — Session 4: Testing, spack setup, libMesh/PETSc install

**Topics:** Running Phase 1 tests, discovering missing HPC libraries, bootstrapping spack.

### Actions taken

**Phase 1 tests (Python):**
- `python3 -m pytest tests/test_mesh.py` — 9/9 pass after fixing `test_no_degenerate_elements`
  (already committed as `1a48a00` in Session 3; confirmed clean on pull)

**Library discovery:**
- libMesh and PETSc are no longer available as Oscar system modules
- No pre-existing user install found under `~/` or anywhere in `/oscar/rt`
- System spack tree at `/oscar/rt/9.6/25/spack` has no spack CLI and no petsc/libmesh packages

**Spack bootstrap:**
- Cloned spack 0.23.1 to `~/data/spack` (more quota headroom than `~/`)
- Compiler registered: gcc@11.5.0, os=oracle9 (`~/.spack/linux/compilers.yaml`)
- Externals registered (`~/.spack/packages.yaml`):
  - `openmpi@4.1.8` — `/oscar/rt/9.6/25/spack/x86_64_v3/openmpi-4.1.8-…`
  - `cmake@3.26.5` — `/usr`
- Activation: `source ~/data/spack/share/spack/setup-env.sh`

**libMesh + PETSc install (in progress):**
- Spec: `libmesh@1.7.1+exodusii+mpi+petsc%gcc@11.5.0`
- PETSc will include: Hypre (BoomerAMG), SuperLU-dist, HDF5, METIS/ParMETIS
- Build launched with `spack install --jobs 8 libmesh+exodusii+mpi+petsc`
- Log: `~/data/spack_install.log`
- Session ended while build was still running

### Next steps
- Confirm spack build completes successfully
- `spack load libmesh` and verify `libmesh-config` is on PATH
- `cmake -DLIBMESH_DIR=$(spack location -i libmesh) -B build .` and build project
- Run `ctest` to execute `test_spaces` (Phase 2 C++ test)

---

## 2026-03-30 — Session 5: spack fresh install, libMesh/PETSc completed, test_spaces passing

**Topics:** Completed PETSc + libMesh install from scratch; debugged concretization and
runtime issues; all C++ tests passing.

### Actions taken

**Spack upgrade and clean slate:**
- Previous session's install was killed on the login node (openblas build OOM)
- Spack had been upgraded to 1.1.1 (separate packages repo `releases/v2025.11`)
- Ran `spack uninstall --all -y` to clear the stale DB; install tree was already empty

**Externals extended (`~/.spack/packages.yaml`):**
- Added `hdf5@1.14.6+mpi` at `/oscar/rt/9.6/25/spack/x86_64_v3/hdf5-1.14.6-…`
  (available as `module load hdf5/1.14.6-663r`)
- Also registered: openblas@0.3.29, netlib-lapack@3.12.1, metis@5.1.0, python@3.12.8

**Spack 1.x concretization bug — workaround:**
- Spec `petsc@3.24.1+hypre` failed: "Cannot satisfy 'hypre@2.32:' (2.12.1)"
- Root cause: petsc package.py has `depends_on("hypre+mpi~internal-superlu", when="+hypre")`;
  the `internal-superlu` variant only exists for `hypre@:2.12.1`; the clingo solver in
  spack 1.x treated `~internal-superlu` as requiring that variant to exist, pinning
  hypre to ≤2.12.1 which conflicts with petsc's `hypre@2.32:` lower bound.
- Fix: patched local packages repo (`~/.spack/package_repos/.../petsc/package.py`),
  replacing `depends_on("hypre+mpi~internal-superlu", when="+hypre")` with
  `depends_on("hypre+mpi", when="+hypre")` (the `~internal-superlu` constraint is
  vacuous for modern hypre and was causing the solver to fail).

**PETSc 3.24.1 installed** (latest in channel):
- Spec: `petsc@3.24.1+hypre+hdf5+metis+superlu-dist+mpi%gcc@11.5.0`
- Built from source: hypre@3.0.0, superlu-dist@9.1.0, parmetis@4.0.3
- System externals reused: openmpi, openblas, hdf5, metis, cmake, python
- Install path: `/users/pcazeaux/data/spack_store/linux-cascadelake/petsc-3.24.1-omxvdtjq…`

**libMesh 1.7.6 installed** (latest in channel):
- Spec: `libmesh@1.7.6+exodusii+mpi+petsc%gcc@11.5.0 ^petsc@3.24.1+…`
- Install path: `/users/pcazeaux/data/spack_store/linux-cascadelake/libmesh-1.7.6-twdoaauu…`

**CMake build clean:**
- `spack load petsc libmesh` → sets `PETSC_DIR` and `CMAKE_PREFIX_PATH`
- `cmake -B build -S .` → finds `libmesh-config`, all includes resolved
- `cmake --build build -j4` → clean build, all targets

**CTest fixes:**
- Tests were launched directly (not via mpirun); libMesh uses MPI and calls
  `MPI_Init` internally → failed with "OPAL ERROR: Unreachable … direct launched using srun"
- Fix: changed `add_test(… COMMAND ${name})` → `add_test(… COMMAND mpirun -n 1 ${name})`
  in `tests/CMakeLists.txt`

**mesh.all_second_order() fix:**
- After mesh fix, test crashed: "Unsupported order: 257" in `DofMap::reinit` /
  `fe_lagrange.C:669` during `EquationSystems::init()`
- Root cause: Gmsh generates TRI3 (linear triangle) elements by default; libMesh's
  DofMap requires TRI6 (6-node quadratic triangle) elements to distribute SECOND-order
  Lagrange DOFs.  The value 257 is an uninitialized/invalid `Order` enum encountered
  when the DofMap iterates over TRI3 elements for a P2 variable.
- Fix: added `mesh.all_second_order()` immediately after `mesh.read()` in both
  `src/main.cpp` and `tests/test_spaces.cpp`; this upgrades all TRI3→TRI6 in-place
  by inserting edge midpoint nodes.

**All tests pass:**
```
1/1 Test #1: test_spaces ........ Passed  0.51 sec
100% tests passed, 0 tests failed out of 1
```

### Files modified this session
- `src/main.cpp` — added `mesh.all_second_order()`
- `tests/test_spaces.cpp` — added `mesh.all_second_order()`
- `tests/CMakeLists.txt` — wrap test command with `mpirun -n 1`
- `notes/conversation_log.md` — this entry

### Next steps (Phase 3)
- Configure PETSc SNES/KSP: FGMRES (restart=100) + BoomerAMG via command-line flags
- Implement Stokes → Navier-Stokes initialization sequence in `main.cpp`
- `tests/test_stokes.cpp` — verify Stokes solve converges, residual < 1e-8

---

## 2026-03-30 — Session 6: Phase 3 solver review and test_stokes

**Topics:** Reviewed Phase 3 solver configuration code (written in a prior sub-session);
discussed parameter choices; wrote and ran `tests/test_stokes.cpp`.

### Solver parameter review

The Phase 3 modifications to `src/main.cpp`, `src/channel_flow_system.h`, and
`src/channel_flow_system.cpp` were inspected and the following design decisions were
confirmed or flagged:

**FGMRES over standard GMRES:**
- Correct choice: BoomerAMG applies a different linear solve on each application,
  making it a variable preconditioner.  Standard GMRES loses orthogonality guarantees
  with variable preconditioners; FGMRES handles this correctly.

**GMRES restart = 100:**
- Reasonable for O(10k–50k) DOF saddle-point systems.  Too small → slow convergence
  from frequent restarts; too large → O(restart²) memory growth.

**KSP_RTOL = 1e-10 (same as SNES_RTOL):**
- This is tight. Newton does not need the linear subproblem solved to full accuracy
  on early iterations (inexact Newton / Eisenstat-Walker strategy would be better).
  However, the current setting is not incorrect; it just does extra work on early
  Newton steps.  Revisit if solve times are too long on finer meshes.

**BoomerAMG on the full saddle-point block:**
- Not the textbook-optimal choice (block Schur / fieldsplit preconditioners are
  better for large Re or fine meshes), but empirically effective at low Re (≤ 20)
  where the system is viscosity-dominated.  Accepted for now; will revisit if
  KSP iteration counts become large.

**Stokes → Navier-Stokes initialization:**
- Solves the linear Stokes problem first (set_stokes_mode(true)) to obtain a
  physically reasonable initial velocity field, then switches to full NS for the
  Newton iteration.  Prevents Newton divergence from a zero initial guess, especially
  at Re > 1.  Stokes is linear so Newton takes exactly one step.

**Pressure null space:**
- With do-nothing (Neumann) outlet and Dirichlet conditions elsewhere, the pressure
  is determined only up to an additive constant.  PETSc/libMesh may or may not
  handle this automatically.  If the Stokes solve produces a pressure that drifts
  or KSP stalls, a pressure pin (fix p = 0 at one outlet node) should be added.
  Flagged as a potential issue to watch during test execution.

### Decisions made

- No changes to existing Phase 3 code; confirmed correct as written.
- `tests/test_stokes.cpp` written to verify: (1) Newton converges in 1 step,
  (2) final nonlinear residual < SNES_ATOL, (3) no-slip BCs satisfied on walls
  and cylinder, (4) inlet u-velocity is non-zero (sanity).
- `tests/CMakeLists.txt` updated: `add_solver_test(test_stokes ...)` uncommented.

### Files created/modified this session
- `tests/test_stokes.cpp` — new
- `tests/CMakeLists.txt` — uncommented test_stokes
- `PLAN.md` — Phase 3 marked done
- `notes/conversation_log.md` — this entry

### Solver debugging attempts during test_stokes development

Three solver configurations were tried before landing on a working test setup.
Recorded here in full detail; the null space issue is revisited in Session 7.

**Attempt 1 — FGMRES + BoomerAMG (production config)**

PETSc options set:
```
-ksp_type fgmres  -ksp_gmres_restart 100
-ksp_rtol 1e-10   -ksp_max_it 500
-pc_type hypre    -pc_hypre_type boomeramg
```
Result: KSP hit the 500-iteration ceiling with residual only reduced from
0.0252 → 0.00436 (factor ≈ 6 over 500 steps — essentially stalled).  Newton
then attempted backtracking but every step was rejected; eventually
"Inexact Newton step FAILED at step 3" → `MPI_ABORT`.

Root-cause diagnosis: the saddle-point system has a pressure null space.  With
do-nothing (Neumann) outlet, pressure is determined only up to a constant.
BoomerAMG receives a singular (or near-singular) matrix and cannot find a good
coarsening hierarchy; the iteration count grows without bound.  The null space
was flagged as a risk in the pre-session review but confirmed here empirically.

**Attempt 2 — `preonly` + LU (direct solve)**

```
-ksp_type preonly
-pc_type lu
```
Result: PETSc error at runtime:
```
KSP of type preonly doesn't make sense with nonzero initial guess
you probably want a KSP of type Richardson
```
Root cause: libMesh calls `KSPSetInitialGuessNonzero(ksp, PETSC_TRUE)` before
every linear solve so it can reuse the previous Newton iterate as a warm start.
`preonly` applies the preconditioner exactly once with no iteration, so it
cannot accept a nonzero initial guess (it would just overwrite it with one
preconditioner application and declare "done", producing a wrong answer).

**Attempt 3 — GMRES + LU (direct solve, correct form)**

```
-ksp_type gmres  -ksp_rtol 1e-12  -ksp_atol 1e-12  -ksp_max_it 1000
-pc_type lu
```
Result: converged in 1 GMRES step (residual 1.5e-14); Newton converged to
9.9e-17 in 1 outer iteration.  Test passes.

Why this works: with LU as the exact preconditioner, the preconditioned system
is the identity (Ax = b → M⁻¹Ax = M⁻¹b with M=A → solve in 1 step).  GMRES
handles the nonzero initial guess correctly (it subtracts the initial residual
and solves the correction equation).  For the 6k-DOF test mesh, LU is cheap.

Bug also found in `max_abs_on_boundary`: `dof_indices(elem, …)` returns all
6 TRI6 DOFs per element (including interior nodes with nonzero velocity);
should query only the 3 DOFs on the boundary side.  Fixed by using
`elem->build_side_ptr(s)` to get the side element, then calling
`dof_map.dof_indices(side.get(), dofs, var)`.

**Status after Session 6 commit:**
The test uses GMRES+LU (robust, tests physics).  Production `main.cpp` still
uses FGMRES+BoomerAMG.  The null space / preconditioner issue is a known gap
that must be resolved before production runs.

### Next steps (Session 7 — null space / preconditioner fix)
- Add single-node pressure pin at outlet corner (BID 5) to eliminate the null space
- Tag node in `main.cpp`, tests after mesh load; add DirichletBC in `init_data()`
- Re-test production config FGMRES+BoomerAMG with pressure pin
- Update `tests/test_stokes.cpp` to use FGMRES+BoomerAMG (match production)

---

## 2026-03-30 — Session 7: Pressure pin + preconditioner investigation

**Topics:** Implementing pressure null-space fix; diagnosing BoomerAMG failure;
settling on ILU(2) preconditioner for the saddle-point system.

### Pressure pin implementation

Added a single-node p=0 Dirichlet condition at the outlet-bottom corner
(x=CHANNEL_LENGTH, y=0) to eliminate the pressure null space.

**Design:**
- `ChannelFlowSystem::tag_pressure_pin(mesh)` — static method that finds the
  local node nearest to the outlet corner and tags it `BID_PRESSURE_PIN = 5`.
  Must be called after `mesh.all_second_order()` and before `es.init()`.
- `init_data()` registers `DirichletBoundary({BID_PRESSURE_PIN}, {_p_var}, zero)`.
- `BID_PRESSURE_PIN = 5` added to `params.h`.
- `tag_pressure_pin` called in `main.cpp`, `tests/test_stokes.cpp`,
  `tests/test_spaces.cpp`.

**Complication discovered:** libMesh's `DofMap::check_dirichlet_bcid_consistency`
enforces that every DirichletBC boundary ID must actually exist in the mesh's
`BoundaryInfo`.  Tests that construct `ChannelFlowSystem` (including `test_spaces`)
must therefore call `tag_pressure_pin` even if they don't invoke `solve()`.
`test_spaces` was failing with "Could not find Dirichlet boundary id 5 in mesh!"
until the call was added.

### BoomerAMG re-diagnosis with pressure pin

Even with the pressure pin active (null space removed), BoomerAMG still stalled:
KSP hit 500 iterations with residual only reduced from 0.025 → 0.003
(factor ≈ 8 over 500 steps — essentially no progress compared to the required
tolerance of ~4e-8).

Root cause (deeper than null space): the Stokes saddle-point matrix has the form
```
[A  B^T]   A = ν∇²  (positive definite)
[B  0  ]   zero diagonal in the pressure block
```
Standard AMG coarsening relies on the diagonal dominance and sign properties of
the matrix rows (or equivalently, assumes M-matrix/symmetric positive definite
structure).  The pressure block rows have *zero diagonal*; AMG cannot build a
valid coarsening hierarchy for them regardless of null-space treatment.

The principled solution is a **fieldsplit preconditioner**: apply AMG only to
the (1,1) velocity block and approximate the pressure Schur complement
S = −B A⁻¹ B^T separately (e.g., with a pressure-Laplacian approximation).
This is deferred as a future improvement for large meshes/high Re.

### Preconditioner decision: ILU(2)

Switched to `ILU(k)` (incomplete LU with fill level k=2) for the full coupled
system.  ILU is a purely algebraic preconditioner: it factorizes the matrix
without assumptions on diagonal dominance, so the saddle-point structure is
handled correctly.

Parameters in `params.h`:
```cpp
constexpr int ILU_FILL = 2;  // pc_factor_levels
```
PETSc options (in `main.cpp` and tests):
```
-ksp_type fgmres  -ksp_gmres_restart 100
-ksp_rtol 1e-10   -ksp_max_it 500
-pc_type ilu      -pc_factor_levels 2
```

**Observed behaviour on coarse mesh (6078 DOFs):**
- KSP: 500 inner iterations, final residual 4.8e-33 (well below tolerance 4e-8)
- Newton: 1 outer step (correct for linear Stokes), final residual 2.7e-17
- Wall time: ≈ 1.4 s including mesh load

500 iterations is higher than ideal for a well-preconditioned system.
Stokes saddle-point conditioning limits ILU(2) efficiency; higher fill (ILU(4))
or a block preconditioner would reduce the count.  For the targeted problem sizes
(≤ 200k DOFs, Re ≤ 20) this is acceptable, but a fieldsplit preconditioner should
be implemented before attempting very fine production meshes.

### Files created/modified this session
- `src/params.h` — added `BID_PRESSURE_PIN = 5`, `ILU_FILL = 2`; updated comment
- `src/channel_flow_system.h` — added `tag_pressure_pin` declaration
- `src/channel_flow_system.cpp` — implemented `tag_pressure_pin`; added
  pressure-pin DirichletBC in `init_data()`; added `#include boundary_info/node`
- `src/main.cpp` — call `tag_pressure_pin`; switched to ILU(2) options
- `tests/test_spaces.cpp` — call `tag_pressure_pin`
- `tests/test_stokes.cpp` — call `tag_pressure_pin`; switched to ILU(2) options
- `notes/conversation_log.md` — this entry

### Next steps (Phase 4)
- ExodusII output: write velocity/pressure fields to `results/channel_flow.e`
- Drag/lift post-processing: boundary integral over cylinder surface (BID 4)
- `tests/test_output.cpp` — verify file written, readable, drag finite, lift ≈ 0 (Stokes)
- Future: fieldsplit preconditioner (AMG on velocity block + Schur for pressure)

---

## 2026-03-30 — Session 8: KSP convergence diagnosis and std::to_string bug

**Topics:** Diagnosing why KSP reported 500 iterations and DIVERGED_ITS despite the
residual appearing to satisfy the tolerance; finding and fixing the root cause.

### The false convergence story (retracted)

After Session 7 the KSP output showed "Linear solve finished, step 500, residual
4.80055e-33" and "DIVERGED_ITS".  I initially offered incorrect explanations:

1. **First claim (wrong):** "ILU(2) needs 500 iterations for a 6k saddle-point
   system; the residual happened to be near zero at the last step because FGMRES
   with restart=100 only checks convergence at restart boundaries."

   The user correctly pointed out: if the residual at step 500 is essentially zero,
   it cannot be above the tolerance at step 499 (restarted FGMRES is monotone
   *within* each cycle).

2. **Second claim (wrong):** "The residual might have been above the tolerance at
   step 499 and fell sharply at step 500 because of how the Krylov subspace fills."

   The user added: the problem is linear (Stokes mode), so we should have exactly
   one Newton iteration; the claim about FGMRES stagnation behavior was speculative
   and not supported.

The user asked to add `-ksp_monitor -ksp_converged_reason`.  The monitor trace
showed the residual satisfying the 3.97e-8 tolerance at step 26 (1.76e-8 < 3.97e-8)
yet the solver continued to step 500.  Combined with `-ksp_view` showing
`tolerances: relative=0., absolute=1e-50`, the actual root cause became clear.

### Root cause: std::to_string truncates small floats to zero

```cpp
std::to_string(1e-10)  // produces "0.000000"  (printf %f, 6 decimal places)
std::to_string(1e-8)   // produces "0.000000"
```

All floating-point tolerances passed to `PetscOptionsSetValue` via `std::to_string`
were silently converted to the string "0.000000", which PETSc read as rtol=0 and
atol=0.  PETSc's convergence criterion was therefore never satisfiable by residual,
and the solver always ran to `max_it` (500), exiting with DIVERGED_ITS.

The integer parameters (`GMRES_RESTART=100`, `KSP_MAX_IT=500`, `ILU_FILL=2`) were
unaffected because `std::to_string(int)` produces the correct decimal string.
This is why `max_it=500` was respected but `ksp_rtol=1e-10` was not.

### Fix

Replaced `std::to_string(double)` with a `petsc_str()` helper in `main.cpp` and
`test_stokes.cpp`:
```cpp
static std::string petsc_str(double v) {
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(6) << v;
    return ss.str();
}
```
`petsc_str(1e-10)` → `"1.000000e-10"`, which PETSc reads correctly.

### Verified behaviour after fix

```
35 KSP Residual norm 1.51e-12
Linear solve converged due to CONVERGED_RTOL iterations 35
```
The Stokes linear system now converges in **35 FGMRES iterations** (down from 500),
with `CONVERGED_RTOL`.  ILU(2) is an effective preconditioner for this problem
once the tolerance is actually communicated to PETSc.

### Files modified this session
- `src/main.cpp` — replaced `std::to_string(KSP_RTOL)` with `petsc_str()`
- `tests/test_stokes.cpp` — same fix; added `petsc_str()` helper

### Next steps (Phase 4)
- ExodusII output and drag/lift coefficients
- Add `petsc_str` to a shared utility header to avoid duplication

---

## 2026-03-30 — Session 9: Fieldsplit AMG Preconditioner

**Topics:** Block Schur complement preconditioner; programmatic PETSc IS configuration;
`petsc_str` consolidation; debugging `pc_fieldsplit_schur_pre_type` options parsing.

### Motivation

ILU(2) converges in ~35 iterations on the 6k-DOF coarse mesh but its iteration
count grows with mesh size.  For production runs (50k–200k DOFs) a scalable
preconditioner is needed.

### Implementation

The standard scalable preconditioner for velocity-pressure saddle-point systems
is the block Schur complement (Silvester & Wathen 1994):

```
P = [A   0]     A → BoomerAMG (mesh-independent V-cycle)
    [B  -S̃]     S̃ → Jacobi on Sp (assembled Schur approximation)
```

PETSc's `-pc_type fieldsplit` supports this, but it requires the velocity and
pressure DOF index sets to be registered programmatically via `PCFieldSplitSetIS`
(cannot be done through `PetscOptionsSetValue` alone).

**New function:** `ChannelFlowSystem::configure_fieldsplit(sys, mesh, newton)`
- Calls `newton.reinit()` to ensure the linear solver is created
- Casts to `PetscLinearSolver`, calls `.ksp()` to get the KSP handle
- Iterates mesh nodes and partitions DOFs into velocity (u+v) and pressure index sets
- Calls `PCFieldSplitSetIS`, `PCFieldSplitSetSchurFactType`, `PCFieldSplitSetSchurPre`
  programmatically (not via options strings, which were silently ignored)
- Sub-PC options (`-fieldsplit_velocity_pc_type hypre`, etc.) still set via
  `PetscOptionsSetValue` and picked up at `PCSetUp` time

**New file:** `src/petsc_utils.h` — moves the `petsc_str()` helper out of
`main.cpp` and `test_stokes.cpp` to eliminate duplication.

### Key debugging findings

1. **`pc_fieldsplit_schur_pre_type a11/selfp` was silently ignored.**
   PETSc's `PCFieldSplitSetFromOptions` only reads `schur_pre_type` when
   `PC_FIELDSPLIT_SCHUR_FACT_*` is already set to a Schur type at the time
   options are processed.  The option was always reported as "unused".
   Fix: call `PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL)`
   directly in C++.

2. **`DIAG` vs `LOWER` factorisation.**
   - `PC_FIELDSPLIT_SCHUR_FACT_DIAG`: treats velocity and pressure blocks
     independently; 176 FGMRES iterations on the coarse mesh.
   - `PC_FIELDSPLIT_SCHUR_FACT_LOWER`: lower-triangular block factorisation,
     captures the velocity-pressure coupling term B*A⁻¹ applied to the residual
     before the pressure correction; **35 iterations** — matches ILU(2) exactly.

3. **ILU for pressure sub-PC broke convergence.**
   The assembled Sp = A10 Diag(A00)⁻¹ A01 has negative diagonal entries (our
   continuity sign convention: -∫q∇·u dx → B = -divergence).  ILU factorisation
   encountered negative pivots and diverged.  Jacobi on the negative diagonal is
   mathematically stable (scales by -1/diagonal).

4. **`-ksp_rtol` override was preventing convergence.**
   libMesh's `NewtonSolver` computes an adaptive linear solve tolerance each step
   (e.g., 3.97e-08 for the Stokes step).  Our `PetscOptionsSetValue("-ksp_rtol",
   "1e-10")` was applied AFTER libMesh's `KSPSetTolerances` call (via
   `KSPSetFromOptions` in `PetscLinearSolver::solve()`), overriding the adaptive
   tolerance with a much tighter one (1e-10), forcing the solver to run to max_it.
   Fix: removed the `-ksp_rtol` option; linear tolerance is now set exclusively
   by libMesh's inexact-Newton framework.  `KSP_RTOL` removed from `params.h`.

### Result

Both tests pass.  `test_stokes`: 35 FGMRES iterations, 1 Newton step,
`CONVERGED_RTOL`.  The preconditioner is mesh-scalable (BoomerAMG + Schur).

### Files created/modified
- `src/petsc_utils.h` — new: `petsc_str()` helper (eliminates duplication)
- `src/channel_flow_system.h` — added `configure_fieldsplit` declaration
- `src/channel_flow_system.cpp` — implemented `configure_fieldsplit` (~60 lines)
- `src/params.h` — removed `ILU_FILL` and `KSP_RTOL`; updated solver comment
- `src/main.cpp` — replaced ILU options with fieldsplit options; added
  `configure_fieldsplit` call; switched to `petsc_utils.h`
- `tests/test_stokes.cpp` — same options update; updated iteration check to ≤ 100

### Next steps (Phase 4)
- ExodusII output (write velocity/pressure fields in `.e` format)
- Drag and lift coefficient computation on the cylinder boundary
- `test_output.cpp`

---

## 2026-03-30 — Session 10: Fix Fieldsplit Sub-PC for Steady Navier-Stokes

**Topics:** BoomerAMG failure diagnosis on Oseen operator; correct AMG target for
incompressible NS; swap velocity/pressure sub-PC choices.

### Problem

After Session 9 implemented fieldsplit, `test_stokes` passed (35 iterations) but
the full NS Newton solve (`Phase 3b`) completely stalled: FGMRES residual made
essentially zero progress in 200 iterations per Newton step.

### Root cause

The velocity block for **steady** NS is the Oseen operator `A = νK + N(u)`.
Without a time-stepping mass matrix term `M/dt`, at Re=20 on the coarse mesh the
cell-Péclet number is ≈ U·h/ν ≈ 0.3·0.1/0.001 = 30.  The block is strongly
non-symmetric and convection-dominated.  Classical BoomerAMG (Falgout coarsening,
symmetric SOR smoother) is designed for SPD systems and breaks down completely here.

The Stokes solve worked because the velocity block was the scalar Laplacian `νK`
(SPD), which AMG handles well.

### Correct preconditioner split

From Elman, Silvester & Wathen "Finite Elements and Fast Iterative Solvers":
AMG is the right tool for the **pressure block** (Sp is SPD), not the velocity
block.  For time-dependent NS the `M/dt` term regularises the velocity block and
makes AMG viable there; for steady NS it is not.

| Block | Preconditioner | Reason |
|-------|---------------|--------|
| Velocity A (Oseen) | ILU(1) | Robust for non-symmetric; not AMG-friendly |
| Pressure Schur Sp | BoomerAMG | SPD assembled matrix; ideal AMG target |

### Result

After swapping the sub-PC options, the full NS Newton solve converges in **4 Newton
steps** on the coarse mesh.  `test_stokes` still passes (61 iterations, 1 Newton
step).

### Files modified
- `src/main.cpp` — velocity: `ilu(1)`, pressure: `hypre boomeramg`
- `tests/test_stokes.cpp` — same
- `src/params.h` — updated solver comment

### Next steps (Phase 4)
- ExodusII output (write `.e` file)
- Drag and lift coefficient computation on cylinder boundary (BID 4)
- `test_output.cpp`

---

## 2026-03-30 — Session 11: Phase 4 — ExodusII Output + Drag/Lift

**Topics:** ExodusII output via libMesh; stress-tensor boundary integration for
drag and lift; sign convention fix; normal convention clarification.

### Implementation

New file `src/drag_lift.cpp` integrates the Cauchy stress traction over the
cylinder boundary (BID_CYLINDER):

```
F = ∫_Γ σ · n_body dS,   σ = -pI + ν(∇u + ∇uᵀ)
```

Iterates over active local elements, detects cylinder-boundary sides via
`has_boundary_id`, calls `fe->reinit(elem, side)` for side FE quadrature,
and manually interpolates u, v, p from the DOF values using phi/dphi.

**Sign convention fix**: libMesh `get_normals()` returns the outward normal of
the fluid element (pointing fluid→cylinder).  Schafer-Turek defines forces
with the outward body normal (cylinder→fluid).  Fix: negate the accumulated
integrals before returning from `compute_drag_lift`.

**ExodusII output**: `libMesh::ExodusII_IO::write_equation_systems` writes the
full velocity/pressure fields in a format readable by ParaView.

### Results on coarse mesh

| Quantity | Value | Schafer-Turek (Re=20) |
|----------|-------|-----------------------|
| C_D (NS)  | 5.52 | 5.57–5.59 |
| C_L (NS)  | 0.011 | 0.010–0.011 |

Good agreement with the benchmark on the coarse mesh (~1% error on C_D).

### Files created/modified
- `src/drag_lift.h` — new: declaration of `compute_drag_lift`
- `src/drag_lift.cpp` — new: boundary stress integration (~115 lines)
- `src/main.cpp` — added ExodusII write + drag/lift printout
- `src/CMakeLists.txt` — added `drag_lift.cpp` to `SOLVER_SOURCES`
- `tests/test_output.cpp` — new: Phase 4 test (file write/read + drag/lift checks)
- `tests/CMakeLists.txt` — enabled `test_output`

### Next steps (Phase 6)
- Python driver (`scripts/run_simulation.py`): CLI arg parsing, SLURM job script
  generation, mesh generation invocation, dry-run mode
- `tests/test_driver.py`

---

## Session 12 — 2026-03-30

### Topics
- Committed Phase 4 work (Session 11 changes)
- Implemented Phase 6: Python driver

### Decisions
- **Driver design**: `scripts/run_simulation.py` with three main functions:
  `parse_args()`, `generate_mesh()`, `generate_slurm_script()`.
  `main()` returns `(args, mesh_path, script)` for testability.
- **Dry-run mode**: `--dry-run` prints SLURM script and mesh command but creates
  no files and does not call `sbatch`.
- **Mesh generation skipped when `--mesh-file` provided**: avoids requiring gmsh
  when the user already has a mesh.
- **SLURM defaults**: 1 node, 1 task, 1h wall time, `batch` partition — conservative
  defaults for the Oscar cluster.

### Files created/modified
- `scripts/run_simulation.py` — new: Python driver (~140 lines)
- `tests/test_driver.py` — new: 10 unit tests covering CLI parsing, SLURM
  structure, mesh invocation, dry-run mode
- `PLAN.md` — marked Phases 4 and 6 as Done

### Test results
- All 10 Phase 6 tests pass (`python -m pytest tests/test_driver.py -v`)

### Next steps (Phase 7)
- Validation against Schafer-Turek benchmark on fine mesh
- `tests/test_validation.py`: C_D ∈ [5.57, 5.59], |C_L| < 0.02 at Re=20

---

## Session 12 (continued) — 2026-03-30

### Topics
- Phase 7: Schafer-Turek benchmark validation on progressively refined meshes

### Mesh convergence study

| Mesh | lc_far/lc_cyl | Elems | DOFs | C_D | C_L | Newton |
|------|---------------|-------|------|-----|-----|--------|
| Coarse | 0.1/0.03 | 244 | 1251 | 5.524 | 0.0106 | 4 |
| Medium | 0.05/0.01 | 1276 | 6087 | 5.524 | 0.0103 | 4 |
| Fine | 0.03/0.005 | 3496 | 16332 | 5.561 | 0.0092 | 4 |
| Refined | 0.02/0.003 | 7990 | 36880 | 5.571 | 0.0104 | 4 |

Schafer-Turek reference: C_D ∈ [5.57, 5.59], C_L ∈ [0.010, 0.011].
The refined mesh (lc_far=0.02, lc_cyl=0.003) hits the benchmark interval.

### Decisions
- **Validation mesh**: lc_far=0.02, lc_cyl=0.003 (~8000 elements, 37k DOFs).
  This is the coarsest mesh that reaches the Schafer-Turek C_D interval with P2/P1.
- **test_validation.py**: Python test that runs the solver binary via mpirun,
  parses C_D/C_L from stdout, and checks against benchmark bounds.  Runs once
  in `setUpModule` and caches output for all 4 checks.
- **Mesh generated on demand**: `_generate_mesh_if_needed()` calls
  `generate_mesh.py` if `channel_refine.msh` does not exist.

### Files created/modified
- `meshes/channel_refine.msh` — refined validation mesh (git-ignored)
- `tests/test_validation.py` — new: 4 tests (Newton convergence, C_D, C_L, ExodusII)
- `PLAN.md` — marked Phase 7 as Done

### Test results
- All 4 Phase 7 tests pass in ~12s (includes solver run)

### Status
All 7 phases complete. The solver is validated against the Schafer-Turek
DFG 2D-1 benchmark at Re=20.

---

## Session 13 — 2026-03-30

### Topics
- Committed Phase 4 (Session 11 leftovers) and Phase 6 (Python driver)
- Phase 7: Schafer-Turek benchmark validation — mesh convergence study
- Planned time-dependent NS extension on new branch

### Mesh convergence study (Phase 7)
Ran the steady NS solver on progressively refined meshes to find the resolution
needed to match the Schafer-Turek C_D ∈ [5.57, 5.59] interval:

| Mesh | lc_far/lc_cyl | Elems | DOFs | C_D | C_L |
|------|---------------|-------|------|-----|-----|
| Coarse | 0.1/0.03 | 244 | 1251 | 5.524 | 0.0106 |
| Medium | 0.05/0.01 | 1276 | 6087 | 5.524 | 0.0103 |
| Fine | 0.03/0.005 | 3496 | 16332 | 5.561 | 0.0092 |
| Refined | 0.02/0.003 | 7990 | 36880 | 5.571 | 0.0104 |

The refined mesh (lc_far=0.02, lc_cyl=0.003) is the coarsest that hits the
benchmark interval.  All 4 validation tests pass in ~12s.

### Decisions
- **New branch `time-dependent-ns`**: created from `main` after all 7 steady
  phases were complete and committed.
- **Time-dependent goal**: visualise recirculation eddies developing from rest
  at Re=5–10.  This was the user's original motivation for the project.
- **Time integrator**: BDF2 (or backward Euler) via libMesh `Euler2Solver`.
  BDF2 is second-order accurate and A-stable.
- **Preconditioner change for time-dependent**: the M/dt mass matrix term
  regularises the velocity block, making it diagonally dominant and nearly SPD.
  BoomerAMG becomes viable for the velocity sub-block (unlike steady where
  ILU was needed due to the non-symmetric Oseen operator at cell-Pe~30).
- **Initial condition**: start from rest (u=v=0) rather than Stokes
  initialisation, so the transient eddy development is visible.
- **Time parameters**: dt=0.025s, T_final=8s, output every 10 steps.
  Convective time L/U_mean ≈ 22s at Re=5; 8s is sufficient to see eddies form.

### Files created/modified
- `CLAUDE.md` — updated problem description for time-dependent NS
- `PLAN.md` — rewritten with 7 new phases for time-dependent extension
- `scripts/run_simulation.py` — new (Phase 6, committed earlier this session)
- `tests/test_driver.py` — new (Phase 6)
- `tests/test_validation.py` — new (Phase 7)

### Next steps
- Phase 1: implement `mass_residual()` in `ChannelFlowSystem`
- Phase 2: time-stepping loop in `main.cpp`
- Phase 3: switch velocity sub-PC to BoomerAMG

---

## Session 14 — 2026-03-30

### Topics
- Implemented time-dependent NS solver (Phases 1–6 of time-dependent plan)

### Key implementation details

**Mass matrix (Phase 1)**: Added `mass_residual()` override to
`ChannelFlowSystem`.  Assembles ∫ρ u̇·w dx for velocity DOFs only (no time
derivative in continuity).  Jacobian: standard mass matrix M_ij = ∫φ_i φ_j dx.

**Time stepping (Phase 2)**: Backward Euler via `libMesh::EulerSolver`
(theta=1.0).  Time loop from t=0 to T_FINAL=8s with dt=0.025s.

**Preconditioner (Phase 3)**: Switched velocity sub-PC from ILU(1) to
BoomerAMG.  With M/dt regularising the velocity block, AMG now works (unlike
the steady case where the non-symmetric Oseen operator required ILU).

**Output (Phases 4–5)**: `ExodusII_IO::write_timestep()` appends snapshots
every OUTPUT_INTERVAL=10 steps to a single `.e` file.  C_D(t) and C_L(t)
printed at each output step.

**Reynolds number**: Changed U_MAX from 0.3 (Re=20) to 0.15 (Re=10) for
the time-dependent eddy development.

### Debugging: time-dependent Dirichlet BCs

The main debugging challenge was getting libMesh's Dirichlet BCs to update
at each time step.

**Problem**: Starting from rest (u=v=0) with a ramp-up inlet BC, the Newton
solver reported "residual 0" at every time step — the solution never evolved.

**Root cause**: libMesh's `DirichletBoundary` constraint values are computed
once during `EquationSystems::init()` (at t=0, where ramp=0) and stored in
`DofMap::_primal_constraint_values`.  These are never automatically updated
for time-dependent BCs.  `enforce_constraints_exactly()` applies the STORED
values, not re-evaluated ones.

**Fix**: Call `DofMap::create_dof_constraints(mesh, t)` at each time step
before `enforce_constraints_exactly(sys)`.  This re-evaluates the Dirichlet
BC function at the new time and updates the stored constraint values.

```cpp
sys.time = t;
sys.get_dof_map().create_dof_constraints(mesh, t);
sys.get_dof_map().enforce_constraints_exactly(sys);
sys.update();
sys.solve();
sys.time_solver->advance_timestep();
```

**Inlet ramp**: Smooth sin² ramp over T_RAMP=2s:
  u(y,t) = 4·U_MAX·y·(H−y)/H² · sin²(πt/(2·T_RAMP))

### Simulation results (coarse mesh, Re=10)

C_D(t) evolution:
  t=0.25: 0.45, t=0.50: 1.75, t=0.75: 3.71, t=1.0: 6.11,
  t=1.5: 10.9, t=2.0: 13.0, t>2.0: 13.03 (steady state)

33 ExodusII snapshots written to results/channel_flow.e.
Solver converges in 1 Newton step per time step during ramp-up.

### Files created/modified
- `src/channel_flow_system.h` — added `mass_residual()` declaration
- `src/channel_flow_assembly.cpp` — added `mass_residual()` implementation
- `src/channel_flow_system.cpp` — time-dependent InletVelocityU with sin² ramp
- `src/main.cpp` — rewritten: EulerSolver time loop, create_dof_constraints
- `src/params.h` — added DT, T_FINAL, T_RAMP, OUTPUT_INTERVAL, THETA; U_MAX→0.15
- `tests/test_timestep.cpp` — new: 20-step test (C_D evolution, multi-step output)
- `tests/CMakeLists.txt` — added test_timestep
- `PLAN.md` — marked Phases 1–6 as Done

### Next steps
- Phase 7: update Python driver for time-dependent parameters
- Run on fine mesh for better visualisation

---

## Session 14 (continued) — 2026-03-30

### Topics
- Phase 7: updated Python driver with solver options
- README.md with problem description, parameter guide, and visualisation instructions

### Changes
- `scripts/run_simulation.py` — added `--re`, `--dt`, `--t-final`,
  `--output-interval` arguments
- `tests/test_driver.py` — added `test_solver_args` test (11 tests total)
- `README.md` — new: problem description, parameter tuning guide (Re, cylinder
  position, time stepping, mesh resolution), visualisation instructions
  (ParaView local, client-server, Open OnDemand), project structure, solver
  details, references
- `PLAN.md` — all 7 phases marked Done

### Status
All time-dependent phases complete.  Solver produces ExodusII time series
showing eddy development from rest at Re=10.

---

## Session 15 — 2026-03-31

### Topics
- Vortex shedding branch (Re=100): time integrator, AMG tuning, Newton Jacobian

### Time integrator: Crank-Nicolson (θ=0.5)

Switched from backward Euler (θ=1, first-order) to EulerSolver with θ=0.5
(Crank-Nicolson, second-order).  This allows dt=0.02 instead of dt=0.005
(4× fewer steps) for equivalent temporal accuracy.  Euler2Solver was also
tested but rejected: it evaluates θ·F(u_new) + (1-θ)·F(u_old) with a
Jacobian that only captures the u_new contribution, giving inherently linear
Newton convergence with factor (1-θ) = 0.5.

### AMG tuning

Experimented with Guermond's BoomerAMG recipe:
- **Strong thresholds** (0.1 velocity, 0.7 pressure): 32% fewer FGMRES
  iterations.  Kept.
- **Chebyshev smoother**: degraded convergence on non-symmetric velocity
  block.  Rejected; default hybrid symmetric GS is better here.
- **Fixed V-cycle count** (Richardson KSP, 2 or 4 V-cycles): too weak as
  sub-PC inside Schur fieldsplit — FGMRES hit 200-iteration limit.
  Rejected; keeping preonly (AMG solves to convergence internally).

ILU(1) for velocity was also tested: 22% fewer FGMRES iterations than
default AMG on the coarse mesh, but strong-threshold AMG was even better.

### Newton Jacobian: restoring quadratic convergence

**The problem**: With Crank-Nicolson (θ=0.5), Newton converged linearly
with factor 0.5, taking 10 iterations per time step instead of the expected
2–3 for quadratic convergence.

**Wrong explanation (initially proposed by Claude)**: "The convergence IS
quadratic, but the quadratic constant C is large at Re=100 due to strong
advection."  The user correctly challenged this: at residual values ~1e-7,
we should be well within the quadratic convergence basin.  A test with very
tight inner FGMRES tolerance (1e-10 relative) confirmed the numbers were
identical — the linear solve accuracy was never the issue.

**First partial fix**: Multiplied Jacobian entries in `element_time_derivative`
by `context.get_elem_solution_derivative()` (= θ).  This is the chain-rule
factor ∂u_θ/∂u = θ, required by libMesh's DiffContext API.  Reduced Newton
from 10 to 6 iterations, but convergence was still not quadratic.

**Root cause (found by reading libMesh source)**: Two bugs in `mass_residual`:

1. **Wrong function for reading u̇**: `c.interior_value()` reads from
   `elem_solution`, which the EulerSolver sets to u_θ (the blended state).
   The correct API is `c.interior_rate()`, which reads from
   `elem_solution_rate` = (u_new − u_old)/dt.

2. **Missing 1/dt factor on mass Jacobian**: The Jacobian of M·u̇ w.r.t. u
   is M·(∂u̇/∂u) = M/dt.  The factor `c.get_elem_solution_rate_derivative()`
   = 1/dt must multiply the mass matrix entries.

Both issues were discovered by reading `euler_solver.C` (cloned from GitHub)
and the libMesh `fem_system_ex1/naviersystem.C` example, which uses
`c.interior_rate()` for mass_residual and asserts `elem_solution_derivative
== 1.0` (i.e., only supports backward Euler).

**Result**: With both fixes, Newton converges quadratically — **2 iterations
per time step** (down from 10 without fixes, 6 with the θ fix only).

### Lesson learned

When Newton convergence is unexpectedly slow, do not accept hand-waving
explanations about "large quadratic constants."  At residual ~1e-7 in a
well-conditioned time-stepping problem, quadratic convergence should be
evident.  The correct approach is: (1) verify the inner solve is not limiting
(tighten KSP tolerance); (2) if unchanged, the Jacobian is wrong — read the
framework source to understand what derivative information the time solver
expects.

### Files modified
- `src/channel_flow_assembly.cpp` — sol_deriv in element_time_derivative;
  interior_rate + rate_deriv in mass_residual; element_constraint unchanged
  (EulerSolver calls it at u_new with derivative=1)
- `src/main.cpp` — EulerSolver θ=0.5, strong thresholds 0.1/0.7
- `src/params.h` — U_MAX=1.5 (Re=100), DT=0.02, THETA=0.5, T_RAMP=1.0
- `tests/test_stokes.cpp`, `tests/test_output.cpp` — apply inlet ramp at
  t=T_RAMP for steady tests

### Next steps
- Run full Re=100 simulation and validate vortex shedding
- Commit conversation log

---

## Session 16 — 2026-03-31

### Topics
- Re=100 simulation runs and preconditioner failure diagnosis
- Movie rendering with pvpython + ffmpeg

### Re=100 simulation: AMG failure and ILU fallback

Ran the solver on the coarse and fine meshes. Both failed at t ≈ 0.84s
(near the end of the inlet ramp, Re_eff ≈ 93):

- FGMRES hit 200 iterations with no convergence
- Newton line search shrank to 1.5e-5 with residual stuck at 5.7e-4
- Root cause: cell Péclet number Pe_h = U·h/ν ≈ 50–140 makes the
  velocity block too non-symmetric for BoomerAMG

**Fix**: switched velocity sub-PC from BoomerAMG to ILU(1). ILU is
algebraically robust for non-symmetric systems.

After the switch, both fine and refined meshes completed successfully:
- Fine mesh (3.5k elem, 16k DOFs): C_D_max = 3.206, C_L_max = 0.972
- Refined mesh (8k elem, 37k DOFs): C_D_max = 3.223, C_L_max = 1.000

Both inside or very close to the Schafer-Turek 2D-2 benchmark interval
(C_D_max ∈ [3.22, 3.24], C_L_max ∈ [0.99, 1.01]).

### Preconditioner summary by Reynolds number

| Re | Velocity sub-PC | Reason |
|----|-----------------|--------|
| 5–10 | BoomerAMG (θ_s=0.1) | M/dt dominates, block nearly SPD |
| 20+ | ILU(1) | Advection dominates, block non-symmetric |

This matches the experience from the steady solver (Session 10): AMG fails
on the non-symmetric Oseen operator at moderate-to-high Re.

### Movie rendering

ParaView's pvpython segfaulted on the compute node (no GPU, no display).
Fix: use Xvfb to create a virtual X11 display:

```bash
Xvfb :99 &; DISPLAY=:99 pvpython scripts/render_movie.py
```

The render_movie.py script:
1. Loads ExodusII file, computes velocity magnitude via Calculator filter
2. Sets up 2D camera view with color bar
3. Renders 81 PNG frames
4. Assembles into MP4 with ffmpeg (libx264, 15 fps, CRF 18)

Output: results/channel_flow.mp4 (559 KB, 5.4s).

### Files created/modified
- `src/main.cpp` — velocity sub-PC switched to ILU(1)
- `scripts/render_movie.py` — new: pvpython rendering + ffmpeg assembly
- `notes/solver_notes.md` — new: comprehensive solver/preconditioner notes
- `PLAN.md` — all 6 phases marked Done

### Status
Vortex shedding branch complete. Simulation validated against
Schafer-Turek 2D-2 benchmark. Movie produced.

---

## Session 17 — 2026-03-31

### Topics
DG Navier-Stokes: formulation choice, libMesh DG attempt, switch to deal.II.

### Decisions and rationale

1. **DG formulation chosen: SIPG + Lax-Friedrichs + equal-order P_k/P_k**
   - SIPG (Symmetric Interior Penalty) for viscous terms — simplest DG diffusion form
   - Lax-Friedrichs upwind flux for advection — smooth Jacobian, robust for all Re
   - Equal-order with pressure-jump stabilization (Shahbazi/Cockburn 2007)
   - Starting at P2/P2, can raise order by changing one parameter

2. **libMesh's FEMSystem does NOT support DG assembly**
   - `DGFEMContext` provides data structures (neighbor shape functions,
     4-block Jacobians) but `FEMSystem::assembly()` never calls
     `neighbor_side_fe_reinit()` — verified experimentally:
     `dg_terms_are_active()` is always false unless manually called
   - Even after manual `neighbor_side_fe_reinit()`, the neighbor residual
     and DG Jacobian blocks are NOT scattered into the global system
   - LU direct solve gives `residual inf` — confirming the matrix is
     structurally singular (only volume terms assembled, no face coupling)
   - **Root cause**: FEMSystem was designed for CG; DGFEMContext is for
     user-level DG but requires a custom assembly loop

3. **Switch to deal.II for DG implementation**
   - deal.II has first-class DG support: `MeshWorker::mesh_loop()` with
     automatic `cell_worker`, `face_worker`, `boundary_worker` callbacks
   - `FEFaceValues` handles neighbor side quadrature naturally
   - PETSc backend via `PETScWrappers` (reuse existing PETSc install)
   - Installing deal.II 9.7.1 via spack with minimal dependencies

### Files created/modified
- `src/dg_face_assembly.cpp` — new: DG face integrals (SIPG, upwind, pressure)
  Written and debugged, but unusable due to FEMSystem limitation
- `src/channel_flow_system.h/.cpp` — modified for DG (L2_LAGRANGE, build_context,
  compute_internal_sides, fieldsplit for element-local DOFs)
- `src/channel_flow_assembly.cpp` — minimal changes (volume integrals identical)
- `src/params.h` — added DG_ORDER, SIGMA_V, GAMMA_P; increased KSP_MAX_IT
- `src/main.cpp` — removed all_second_order, tag_pressure_pin, constraint updates
- `tests/test_dg_spaces.cpp` — new: DG space verification (PASSES)
- `tests/test_dg_stokes.cpp` — new: DG Stokes solve (FAILS: singular matrix)
- `PLAN.md` — updated for deal.II approach

### Next steps
1. Wait for deal.II spack install to complete
2. Set up deal.II build infrastructure (CMakeLists.txt)
3. Implement DG solver using deal.II's MeshWorker::mesh_loop()
4. Validate against Schafer-Turek benchmarks

---

## Session 18 — 2026-03-31

### Topics
deal.II installation, IMEX time integration, Schur complement solver design,
Phase 1 implementation (mesh + DOFs + sparsity + mass inverses).

### Decisions and rationale

1. **deal.II 9.7.1 installed via spack**
   - Spec: `+petsc +mpi -threads` (TBB incompatible with GCC 11.5)
   - First build failed on TBB hash<pair> error; disabled threads to fix
   - 43 minutes to build

2. **IMEX time integration (explicit advection, implicit Stokes)**
   - Velocity block A = M_V/dt + θνK_DG becomes SPD (no advection)
   - AMG is optimal for SPD operators → no ILU fallback needed
   - A is constant across all stages and time steps → AMG setup once
   - No Newton iterations needed (system is linear at each stage)
   - IMEX-ARK3 (Kennedy & Carpenter 2003): 3rd order, 4 stages, CFL ≈ 1.7
   - dt ≈ 0.011 s (vs 0.02 fully implicit); 4× more steps but each ~6× cheaper

3. **Pressure Schur complement solver with Cahouet-Chabard preconditioner**
   - Explicitly solve S = C - B^T A^{-1} B via nested Krylov
   - Preconditioner: P_CC = ν·M_p^{-1} + (1/dt)·L_p^{-1}
   - DG advantage: M_V^{-1} and M_p^{-1} are exact (block-diagonal)
   - L_p = B^T M_V^{-1} B assembled explicitly (cheap with DG mass inverse)
   - L_p^{-1} via 2 AMG V-cycles; A^{-1} via AMG (both SPD)

4. **Quad mesh via GridGenerator::channel_with_cylinder**
   - deal.II's FE_DGQ requires quads; our Gmsh mesh was triangles
   - Built-in generator produces exact Schafer-Turek geometry with manifolds
   - Boundary IDs: 0=inlet, 1=outlet, 2=cylinder, 3=walls
   - MappingQ(2) for curved cells near cylinder

5. **Separate DoFHandlers for velocity and pressure**
   - Enables clean block assembly (A, B, C as separate matrices)
   - Velocity: 2 scalar DoFHandlers (u, v share same FE_DGQ)
   - Natural for Schur complement decomposition

### Files created/modified
- `src/dealii/dg_params.h` — centralized parameters (deal.II boundary IDs)
- `src/dealii/dg_navier_stokes.h` — solver class declaration
- `src/dealii/dg_navier_stokes.cc` — Phase 1: mesh, DOFs, sparsity, mass inverses
- `src/dealii/main.cc` — entry point
- `src/dealii/CMakeLists.txt` — build with deal.II (cmake + deal_ii_setup_target)
- `PLAN.md` — updated with IMEX-ARK3 + Schur complement approach
- `notes/solver_notes.md` — detailed IMEX and Cahouet-Chabard analysis
- `notes/conversation_log.md` — this entry

### Phase 1 results
- 432 cells (1 global refinement), 11,664 total DOFs
- 9 DOFs/cell/component (Q2 on quads)
- 432 element-local mass inverse blocks (9×9) computed
- VTK output verified

### Phase 2 results
All block matrices assembled (4.2s total on 432-cell mesh):
- A = M_V + γ·dt·νK_DG: SIPG via MeshWorker::mesh_loop (cell + face + boundary)
- B, B^T: pressure gradient/divergence (volume terms)
- C: pressure-jump stabilization (face integrals via FEInterfaceValues)
- L_p = B^T M_V^{-1} B: element-local assembly using pre-computed M_V^{-1}
- M_V, M_p: mass matrices

Key fixes during Phase 2:
- SparsityPattern for B^T must be stored as class member (ObserverPointer)
- boundary_worker for A uses direct matrix.add() (not CopyData) to avoid
  conflicting with cell_worker CopyData
- FEInterfaceValues handles jump/average operators for SIPG automatically

### Schur complement solver (Phase 5, partial)
- SchurOperator: matrix-free S = C - B^T A^{-1} B via nested CG+ILU
- CahouetChabard: P_CC = ν M_p^{-1} + (1/dt_eff) L_p^{-1}
  - M_p^{-1}: exact element-local (DG block-diagonal)
  - L_p^{-1}: ILU on DG SIPG Laplacian (replaces volume-only B^T M_V^{-1} B
    which was block-diagonal and singular)
- Pre-solve converges in 1 CG iteration (velocity RHS very small)
- GMRES on Schur complement diverges — needs debugging:
  - RHS formulation may be incomplete (missing mass term contribution)
  - L_p boundary conditions (Dirichlet penalty on non-outlet) may be wrong
    for the pressure preconditioner
  - Schur complement operator signs need verification

### Solver debugging (continued in session)

Systematic audit identified 4 bugs:
1. ✅ B only had x-component — fixed: B is now (2·n_scalar) × n_pres
2. ✅ B missing face integrals — fixed: added ∫_F {p}[[w·n]] for interior + boundary
3. S null space — pinned DOF 0 in Schur operator (identity row)
4. ✅ Sign: negated S so -S is PSD, matching the PD preconditioner

Convergence status after fixes:
- Unpreconditioned FGMRES: 0.114 → 0.022 in 100 iters (slowly converging)
- With Cahouet-Chabard: 0.114 → 0.108 in 100 iters (preconditioner ineffective)
- The operator is correct (unpreconditioned converges), but the preconditioner
  is not spectrally equivalent to the Schur complement

Root cause analysis for poor preconditioning:
- dt_eff = γ·dt = 0.004359 → dt_inv = 229, so L_p^{-1} dominates P_CC
- The L_p assembled is a SIPG Laplacian, which may not be spectrally
  equivalent to the actual -S = B^T A^{-1} B - C
- The inner CG solve for L_p^{-1} may need more iterations or better preconditioner
- Need to verify eigenvalue distribution of P_CC^{-1} (-S)

### Bugs found and fixed (continued)
1. **Face accumulation**: `face_matrices.resize(1)` → `emplace_back()`.
   Only last interior face per cell was assembled. Root cause of A singularity.
2. **Double ν in penalty**: `σν²/h` → `σν/h`.
3. **RHS sign**: `+∂w/∂n · u_bc` → `-∂w/∂n · u_bc`.
4. **Spurious pressure RHS**: strong-form B has no boundary continuity terms.
5. **Q_k/Q_k checkerboard**: interior pressure DOFs invisible to div and face
   stabilization. Fix: Q2/Q1 (DG Taylor-Hood), natural inf-sup.

### Steady Stokes verified
- Q2/Q1 DG with UMFPACK: smooth velocity and pressure
- Convergence: per-DOF |u| stable at ~0.21 across 3 refinements
- Strong-form B, no pressure stabilization needed

### Iterative Schur solver (completed)
- Replaced UMFPACK with iterative Schur complement solve
- Velocity: CG + BoomerAMG (strong threshold 0.1, Guermond)
  - Converted deal.II SparseMatrix to PETSc for BoomerAMG
  - 7 CG iterations per component (mesh-independent with AMG)
- Pressure: FGMRES + block-diagonal Cahouet-Chabard
  - P_CC = ν M_p^{-1} + 1/(γdt) L_p^{-1}, both element-local (DG)
  - No AMG needed for pressure preconditioner (DG advantage)
  - 14-22 FGMRES iterations
- Performance: 1.8s total for 38k DOFs in Release mode

### Time stepping (completed, Stokes only)
- IMEX-Euler: explicit advection + implicit Stokes
- Time loop with proper RHS: M_V u_n - dt N(u_n) + coeff × BC terms
- PVD/VTU time series output (single .pvd file for ParaView)
- Velocity and pressure in same VTU (pressure L2-projected to vel mesh)
- Stokes (no advection): stable, smooth, physically reasonable

### Advection operator (IN PROGRESS — blowing up)
- Implemented Lax-Friedrichs DG flux on interior + boundary faces
- BUG: used strong-form volume term `∫(u·∇)u · w` instead of
  conservative DG form `-∫(u⊗u):∇w + ∫ F̂·n · w`
- Strong form doesn't account for discontinuous fluxes at faces
- Causes exponential blowup at any dt (even dt=0.001)
- Fix: rewrite in conservative form with integration by parts

### Next steps
1. Fix advection: conservative DG form with IBP
2. Verify stability with CFL-limited dt
3. Upgrade to IMEX-ARK3 (3rd order, larger CFL)
4. Longer runs, drag/lift computation
5. Validation against Schafer-Turek
