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

    // ── Solve ──────────────────────────────────────────────────────────────
    // Phase 3 will configure KSP/PC options via PETSc command-line flags.
    // Phase 4 will write ExodusII output and compute drag/lift coefficients.
    sys.solve();

    libMesh::out << "Solve complete.\n";
    return 0;
}
