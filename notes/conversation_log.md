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
