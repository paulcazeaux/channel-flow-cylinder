/**
 * @file main.cpp
 * @brief Entry point for the channel-flow solver.
 *
 * Initialises libMesh, loads the mesh, constructs ChannelFlowSystem with a
 * SteadySolver + NewtonSolver chain, and calls solve().
 *
 * Usage:
 *   channel_flow --mesh <path/to/mesh.msh> [--help]
 *
 * Phase 2: FEMSystem with Taylor-Hood P2/P1, Dirichlet BCs applied.
 * Phase 3: PETSc SNES/KSP options set via command line or options file.
 * Phase 4: ExodusII output and drag/lift post-processing.
 */

#include "channel_flow_system.h"
#include "params.h"

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"

#include "petsc_utils.h"

#include <petscsys.h>

#include <iostream>
#include <string>
#include <cstring>

/// @brief Parse --mesh and --help flags from argv; return mesh path or "".
static std::string parse_args(int argc, char** argv)
{
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--help") == 0) {
            std::cout
                << "Usage: channel_flow --mesh <path.msh>\n"
                << "\n"
                << "Solves 2D steady incompressible Navier-Stokes in a\n"
                << "channel past a cylinder (Schafer-Turek 2D-1, Re="
                << Params::RE << ").\n"
                << "\n"
                << "Options:\n"
                << "  --mesh <file>  Path to Gmsh .msh file\n"
                << "  --help         Show this message\n";
            return "";
        }
        if (std::strcmp(argv[i], "--mesh") == 0 && i + 1 < argc)
            return argv[++i];
    }
    // Default to the fine production mesh.
    return "meshes/channel.msh";
}

int main(int argc, char** argv)
{
    // LibMeshInit must be constructed before anything else.
    libMesh::LibMeshInit init(argc, argv);

    const std::string mesh_path = parse_args(argc, argv);
    if (mesh_path.empty())
        return 0;   // --help was printed

    libMesh::out
        << "ChannelFlow: 2D Navier-Stokes past a cylinder\n"
        << "  Re = " << Params::RE << "  U_MAX = " << Params::U_MAX << " m/s\n"
        << "  mesh: " << mesh_path << "\n";

    // ── Load mesh ──────────────────────────────────────────────────────────
    libMesh::Mesh mesh(init.comm());
    mesh.read(mesh_path);
    mesh.all_second_order(); // upgrade TRI3→TRI6 for P2/P1 Taylor-Hood
    ChannelFlowSystem::tag_pressure_pin(mesh); // fix pressure null space
    libMesh::out << "  Elements: " << mesh.n_elem()
                 << "  Nodes: " << mesh.n_nodes() << "\n";

    // ── Build system ───────────────────────────────────────────────────────
    libMesh::EquationSystems es(mesh);
    ChannelFlowSystem& sys = es.add_system<ChannelFlowSystem>("ChannelFlow");

    // Time solver must be set before init(); SteadySolver for steady problems.
    sys.time_solver = std::make_unique<libMesh::SteadySolver>(sys);

    // Nonlinear solver: Newton with backtracking (Phase 3 tunes PETSc flags).
    auto& newton = sys.time_solver->diff_solver();
    newton = std::make_unique<libMesh::NewtonSolver>(sys);
    auto& ns = static_cast<libMesh::NewtonSolver&>(*newton);
    ns.max_nonlinear_iterations  = Params::SNES_MAX_IT;
    ns.absolute_residual_tolerance = Params::SNES_ATOL;
    ns.relative_residual_tolerance = Params::SNES_RTOL;
    ns.verbose = true;

    es.init();

    libMesh::out << "  DOFs: " << sys.n_dofs() << "\n";

    // ── Phase 3: PETSc KSP/PC options ─────────────────────────────────────
    // Outer FGMRES + fieldsplit Schur preconditioner (Silvester & Wathen 1994).
    // Velocity block → BoomerAMG; pressure block → Jacobi on the Schur approx.
    // ISes are registered programmatically via configure_fieldsplit (below).
    PetscOptionsSetValue(NULL, "-ksp_type",           "fgmres");
    PetscOptionsSetValue(NULL, "-ksp_gmres_restart",
                         std::to_string(Params::GMRES_RESTART).c_str());
    // ksp_rtol is left to libMesh's inexact-Newton framework, which sets it
    // adaptively each step.  Overriding with a tight tolerance (1e-10) forces
    // the linear solver to run to the iteration limit.
    PetscOptionsSetValue(NULL, "-ksp_max_it",
                         std::to_string(Params::KSP_MAX_IT).c_str());

    // Fieldsplit type set programmatically in configure_fieldsplit (PCSetType,
    // PCFieldSplitSetType, PCFieldSplitSetSchurFactType, PCFieldSplitSetSchurPre).
    // Sub-PC options are still set here via options string and picked up at PCSetUp.

    // Velocity sub-PC: single BoomerAMG V-cycle
    PetscOptionsSetValue(NULL, "-fieldsplit_velocity_ksp_type",       "preonly");
    PetscOptionsSetValue(NULL, "-fieldsplit_velocity_pc_type",        "hypre");
    PetscOptionsSetValue(NULL, "-fieldsplit_velocity_pc_hypre_type",  "boomeramg");

    // Pressure sub-PC: Jacobi on Sp (the assembled Schur approximation, 706×706)
    PetscOptionsSetValue(NULL, "-fieldsplit_pressure_ksp_type",   "preonly");
    PetscOptionsSetValue(NULL, "-fieldsplit_pressure_pc_type",    "jacobi");

    // Register velocity/pressure IS objects with PETSc fieldsplit.
    ChannelFlowSystem::configure_fieldsplit(sys, mesh, ns);

    // ── Phase 3: Stokes initialisation → Navier-Stokes Newton ─────────────
    // A linear Stokes solve provides a good initial guess and avoids Newton
    // divergence from a zero velocity field (critical at Re > 1).
    libMesh::out << "Phase 3a: Stokes initialisation (linear solve)...\n";
    sys.set_stokes_mode(true);
    sys.solve();

    libMesh::out << "Phase 3b: Navier-Stokes Newton iteration...\n";
    sys.set_stokes_mode(false);
    sys.solve();

    // Phase 4 will write ExodusII output and compute drag/lift coefficients.
    libMesh::out << "Solve complete.\n";
    return 0;
}
