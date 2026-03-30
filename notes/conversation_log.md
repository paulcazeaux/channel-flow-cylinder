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
