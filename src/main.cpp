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
    // FGMRES (restart=100) with ILU(ILU_FILL) on the full saddle-point system.
    // BoomerAMG on the full coupled block does NOT converge: the pressure rows
    // have a zero diagonal, violating AMG coarsening assumptions regardless of
    // null-space handling.  The principled scalable alternative is a fieldsplit
    // preconditioner (AMG on the (1,1) velocity block + Schur complement
    // approximation for pressure), but ILU(k) is robust for the problem sizes
    // targeted here (Re ≤ 20, moderate mesh).
    PetscOptionsSetValue(NULL, "-ksp_type",         "fgmres");
    PetscOptionsSetValue(NULL, "-ksp_gmres_restart",
                         std::to_string(Params::GMRES_RESTART).c_str());
    PetscOptionsSetValue(NULL, "-ksp_rtol",
                         std::to_string(Params::KSP_RTOL).c_str());
    PetscOptionsSetValue(NULL, "-ksp_max_it",
                         std::to_string(Params::KSP_MAX_IT).c_str());
    PetscOptionsSetValue(NULL, "-pc_type",          "ilu");
    PetscOptionsSetValue(NULL, "-pc_factor_levels",
                         std::to_string(Params::ILU_FILL).c_str());

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
