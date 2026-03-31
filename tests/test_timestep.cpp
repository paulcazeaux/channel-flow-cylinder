/**
 * @file test_timestep.cpp
 * @brief Phase 6 test: time-dependent Navier-Stokes with inlet ramp.
 *
 * Pass criteria:
 *   1. Newton converges at each time step (no divergence).
 *   2. Solution evolves between steps (C_D changes during ramp).
 *   3. C_D at the final step is finite and positive.
 *   4. ExodusII output file has multiple time steps.
 *
 * Runs 20 time steps (dt=0.025, covering t=0..0.5) on the coarse mesh.
 */

#include "channel_flow_system.h"
#include "drag_lift.h"
#include "params.h"
#include "petsc_utils.h"

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/euler_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/dof_map.h"

#include <petscsys.h>

#include <sys/stat.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    libMesh::LibMeshInit init(argc, argv);

    const char* mesh_env = std::getenv("TEST_MESH");
    if (!mesh_env) {
        std::cerr << "[test_timestep] FAIL: TEST_MESH not set.\n";
        return 1;
    }

    libMesh::Mesh mesh(init.comm());
    mesh.read(mesh_env);
    mesh.all_second_order();
    ChannelFlowSystem::tag_pressure_pin(mesh);

    libMesh::EquationSystems es(mesh);
    ChannelFlowSystem& sys = es.add_system<ChannelFlowSystem>("ChannelFlow");

    sys.time_solver = std::make_unique<libMesh::EulerSolver>(sys);
    auto& euler = static_cast<libMesh::EulerSolver&>(*sys.time_solver);
    euler.theta = Params::THETA;

    auto& diff_solver = sys.time_solver->diff_solver();
    diff_solver = std::make_unique<libMesh::NewtonSolver>(sys);
    auto& ns = static_cast<libMesh::NewtonSolver&>(*diff_solver);
    ns.max_nonlinear_iterations    = Params::SNES_MAX_IT;
    ns.absolute_residual_tolerance = Params::SNES_ATOL;
    ns.relative_residual_tolerance = Params::SNES_RTOL;
    ns.verbose = false;

    es.init();
    sys.deltat = Params::DT;

    // PETSc options: BoomerAMG for both velocity and pressure blocks.
    PetscOptionsSetValue(NULL, "-ksp_type",           "fgmres");
    PetscOptionsSetValue(NULL, "-ksp_gmres_restart",
                         std::to_string(Params::GMRES_RESTART).c_str());
    PetscOptionsSetValue(NULL, "-ksp_max_it",
                         std::to_string(Params::KSP_MAX_IT).c_str());

    PetscOptionsSetValue(NULL, "-fieldsplit_velocity_ksp_type",         "preonly");
    PetscOptionsSetValue(NULL, "-fieldsplit_velocity_pc_type",          "hypre");
    PetscOptionsSetValue(NULL, "-fieldsplit_velocity_pc_hypre_type",    "boomeramg");

    PetscOptionsSetValue(NULL, "-fieldsplit_pressure_ksp_type",         "preonly");
    PetscOptionsSetValue(NULL, "-fieldsplit_pressure_pc_type",          "hypre");
    PetscOptionsSetValue(NULL, "-fieldsplit_pressure_pc_hypre_type",    "boomeramg");

    ChannelFlowSystem::configure_fieldsplit(sys, mesh, ns);

    // ── Run 20 time steps ─────────────────────────────────────────────────
    constexpr int N_STEPS = 20;
    constexpr double dt = Params::DT;

    (void)mkdir("results", 0755);
    const std::string out_path = "results/test_timestep.e";
    libMesh::ExodusII_IO exo(mesh);

    int output_step = 1;
    exo.write_timestep(out_path, es, output_step, 0.0);

    double cd_first = 0.0, cd_last = 0.0;

    for (int step = 1; step <= N_STEPS; ++step) {
        const double t = step * dt;
        sys.time = t;
        sys.get_dof_map().create_dof_constraints(mesh, t);
        sys.get_dof_map().enforce_constraints_exactly(sys);
        sys.update();

        sys.solve();
        sys.time_solver->advance_timestep();

        // Check 1: Newton converged (no exception thrown).

        // Record drag at first and last output steps.
        if (step == 5 || step == N_STEPS) {
            const auto drag_lift = compute_drag_lift(sys, mesh);
            const double C_D = Params::DRAG_LIFT_NORM * drag_lift.first;
            if (step == 5)  cd_first = C_D;
            if (step == N_STEPS) cd_last = C_D;

            ++output_step;
            exo.write_timestep(out_path, es, output_step, t);
        }
    }

    libMesh::out << "[test_timestep] C_D at step 5:  " << cd_first << "\n";
    libMesh::out << "[test_timestep] C_D at step 20: " << cd_last << "\n";

    // Check 2: Solution evolves — C_D changes during the ramp.
    if (std::abs(cd_last - cd_first) < 1e-6) {
        std::cerr << "[test_timestep] FAIL: C_D did not change between steps.\n";
        return 1;
    }
    libMesh::out << "[test_timestep] C_D evolution check: PASSED.\n";

    // Check 3: C_D at final step is finite and positive.
    if (!std::isfinite(cd_last) || cd_last <= 0.0) {
        std::cerr << "[test_timestep] FAIL: C_D at final step = "
                  << cd_last << " (expected finite, positive).\n";
        return 1;
    }
    libMesh::out << "[test_timestep] C_D finite/positive check: PASSED.\n";

    // Check 4: ExodusII output has multiple time steps.
    struct stat st{};
    if (stat(out_path.c_str(), &st) != 0 || st.st_size == 0) {
        std::cerr << "[test_timestep] FAIL: output file missing or empty.\n";
        return 1;
    }
    if (output_step < 3) {
        std::cerr << "[test_timestep] FAIL: only " << output_step
                  << " output steps (expected >= 3).\n";
        return 1;
    }
    libMesh::out << "[test_timestep] Multi-step output check: PASSED ("
                 << output_step << " steps).\n";

    libMesh::out << "[test_timestep] ALL CHECKS PASSED.\n";
    return 0;
}
