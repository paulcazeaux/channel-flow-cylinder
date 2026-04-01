/**
 * @file main.cpp
 * @brief Entry point for the time-dependent DG channel-flow solver.
 *
 * Initialises libMesh, loads the mesh, constructs ChannelFlowSystem with
 * DG equal-order elements, EulerSolver (Crank-Nicolson) + NewtonSolver,
 * and runs a time-stepping loop from t=0 to T_FINAL.
 *
 * Usage:
 *   channel_flow --mesh <path/to/mesh.msh> [--help]
 *
 * Output: ExodusII time series in results/channel_flow.e, with snapshots
 * every OUTPUT_INTERVAL steps.  Drag/lift coefficients C_D(t), C_L(t)
 * printed at each output step.
 */

#include "channel_flow_system.h"
#include "drag_lift.h"
#include "params.h"

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/euler_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/dof_map.h"

#include "petsc_utils.h"

#include <petscsys.h>

#include <sys/stat.h>
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
                << "Solves 2D time-dependent incompressible Navier-Stokes in a\n"
                << "channel past a cylinder (Schafer-Turek, Re="
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
    return "meshes/channel.msh";
}

int main(int argc, char** argv)
{
    libMesh::LibMeshInit init(argc, argv);

    const std::string mesh_path = parse_args(argc, argv);
    if (mesh_path.empty())
        return 0;

    libMesh::out
        << "ChannelFlow: 2D time-dependent Navier-Stokes past a cylinder\n"
        << "  Re = " << Params::RE << "  U_MAX = " << Params::U_MAX << " m/s\n"
        << "  dt = " << Params::DT << "  T_final = " << Params::T_FINAL << "\n"
        << "  mesh: " << mesh_path << "\n";

    // ── Load mesh ──────────────────────────────────────────────────────────
    libMesh::Mesh mesh(init.comm());
    mesh.read(mesh_path);
    // DG uses L2_LAGRANGE on TRI3 elements — no need for all_second_order().
    libMesh::out << "  Elements: " << mesh.n_elem()
                 << "  Nodes: " << mesh.n_nodes() << "\n";

    // ── Build system ───────────────────────────────────────────────────────
    libMesh::EquationSystems es(mesh);
    ChannelFlowSystem& sys = es.add_system<ChannelFlowSystem>("ChannelFlow");

    // Time solver: EulerSolver with theta=0.5 (Crank-Nicolson).
    // Second-order accurate, A-stable.  Evaluates F at the blended state
    // theta*u_new + (1-theta)*u_old, so Newton sees the full Jacobian
    // and converges quadratically (unlike Euler2Solver which gives linear).
    sys.time_solver = std::make_unique<libMesh::EulerSolver>(sys);
    auto& euler = static_cast<libMesh::EulerSolver&>(*sys.time_solver);
    euler.theta = Params::THETA;

    // Nonlinear solver: Newton with backtracking.
    auto& newton = sys.time_solver->diff_solver();
    newton = std::make_unique<libMesh::NewtonSolver>(sys);
    auto& ns = static_cast<libMesh::NewtonSolver&>(*newton);
    ns.max_nonlinear_iterations  = Params::SNES_MAX_IT;
    ns.absolute_residual_tolerance = Params::SNES_ATOL;
    ns.relative_residual_tolerance = Params::SNES_RTOL;
    ns.verbose = true;

    es.init();

    // Set initial time step for the time solver.
    sys.deltat = Params::DT;

    libMesh::out << "  DOFs: " << sys.n_dofs() << "\n";

    // ── PETSc KSP/PC options ──────────────────────────────────────────────
    PetscOptionsSetValue(NULL, "-ksp_type",           "fgmres");
    PetscOptionsSetValue(NULL, "-ksp_gmres_restart",
                         std::to_string(Params::GMRES_RESTART).c_str());
    PetscOptionsSetValue(NULL, "-ksp_max_it",
                         std::to_string(Params::KSP_MAX_IT).c_str());

    // ── Velocity sub-PC: ASM with ILU(2) ────────────────────────────────
    // DG stiffness has wider stencil from face coupling; use additive
    // Schwarz with ILU(2) on each subdomain.
    PetscOptionsSetValue(NULL, "-fieldsplit_velocity_ksp_type",         "preonly");
    PetscOptionsSetValue(NULL, "-fieldsplit_velocity_pc_type",          "asm");
    PetscOptionsSetValue(NULL, "-fieldsplit_velocity_sub_pc_type",      "ilu");
    PetscOptionsSetValue(NULL, "-fieldsplit_velocity_sub_pc_factor_levels", "2");

    // ── Pressure sub-PC: BoomerAMG ───────────────────────────────────────
    // Laplacian-like Schur complement Sp: strong threshold 0.7.
    PetscOptionsSetValue(NULL, "-fieldsplit_pressure_ksp_type",      "preonly");
    PetscOptionsSetValue(NULL, "-fieldsplit_pressure_pc_type",       "hypre");
    PetscOptionsSetValue(NULL, "-fieldsplit_pressure_pc_hypre_type", "boomeramg");
    PetscOptionsSetValue(NULL, "-fieldsplit_pressure_pc_hypre_boomeramg_strong_threshold", "0.7");

    ChannelFlowSystem::configure_fieldsplit(sys, mesh, ns);
    ChannelFlowSystem::setup_pressure_null_space(sys, mesh, ns);

    // ── Prepare output ────────────────────────────────────────────────────
    (void)mkdir("results", 0755);
    libMesh::ExodusII_IO exo(mesh);

    // Write initial condition (t=0, zero velocity).
    int output_step = 1;
    exo.write_timestep(Params::OUTPUT_FILE, es, output_step, 0.0);
    libMesh::out << "  t = 0  (initial condition written)\n";

    // ── Time loop ─────────────────────────────────────────────────────────
    int step = 0;
    for (double t = Params::DT; t <= Params::T_FINAL + 0.5 * Params::DT;
         t += Params::DT, ++step)
    {
        // Set system time so weak BCs evaluate the ramp at t_new.
        sys.time = t;

        libMesh::out << "\n=== Time step " << step + 1
                     << "  t = " << t << " ===\n";
        sys.solve();
        sys.time_solver->advance_timestep();

        // Output snapshot every OUTPUT_INTERVAL steps.
        if ((step + 1) % Params::OUTPUT_INTERVAL == 0 ||
            t >= Params::T_FINAL - 0.5 * Params::DT)
        {
            ++output_step;
            exo.write_timestep(Params::OUTPUT_FILE, es, output_step, t);

            const auto drag_lift = compute_drag_lift(sys, mesh);
            const double C_D = Params::DRAG_LIFT_NORM * drag_lift.first;
            const double C_L = Params::DRAG_LIFT_NORM * drag_lift.second;
            libMesh::out << "  Output step " << output_step
                         << "  C_D = " << C_D << "  C_L = " << C_L << "\n";
        }
    }

    libMesh::out << "Simulation complete. " << output_step
                 << " snapshots written to " << Params::OUTPUT_FILE << "\n";
    return 0;
}
