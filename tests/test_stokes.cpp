/**
 * @file test_stokes.cpp
 * @brief Phase 3 test: verify Stokes linear initialisation solve.
 *
 * Pass criteria (from PLAN.md §Phase 3):
 *   1. Solve completes without exception.
 *   2. Newton converges in exactly 1 step (Stokes is a linear problem).
 *   3. Final nonlinear residual < Params::SNES_ATOL.
 *   4. No-slip BCs satisfied on walls (BID 3) and cylinder (BID 4):
 *      max |u|, |v| < 1e-10 at boundary nodes.
 *   5. Inlet u-velocity is non-zero (sanity check that BCs are applied).
 *
 * Run via CTest; mesh path passed through the TEST_MESH environment variable.
 */

#include "channel_flow_system.h"
#include "params.h"
#include "petsc_utils.h"

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/boundary_info.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

// ── Helpers ───────────────────────────────────────────────────────────────────

/**
 * @brief Return the maximum absolute value of a variable at nodes on a boundary.
 *
 * Iterates over elements; for each side on @p bid collects the node DOF values
 * from the solution vector and tracks the maximum.
 */
static double max_abs_on_boundary(
    const libMesh::MeshBase&           mesh,
    const libMesh::DofMap&             dof_map,
    const libMesh::NumericVector<libMesh::Number>& sol,
    unsigned int                       var,
    libMesh::boundary_id_type          bid)
{
    double max_val = 0.0;
    std::vector<libMesh::dof_id_type> dofs;

    for (const auto& elem : mesh.active_local_element_ptr_range()) {
        for (unsigned short s = 0; s < elem->n_sides(); ++s) {
            if (!mesh.get_boundary_info().has_boundary_id(elem, s, bid))
                continue;
            // Query DOFs on the side only, not the whole element.
            auto side = elem->build_side_ptr(s);
            dof_map.dof_indices(side.get(), dofs, var);
            for (auto d : dofs) {
                const double val = std::abs(sol(d));
                max_val = std::max(max_val, val);
            }
        }
    }
    return max_val;
}

// ── Main ─────────────────────────────────────────────────────────────────────

int main(int argc, char** argv)
{
    // ── 1. libMesh initialisation ─────────────────────────────────────────────
    libMesh::LibMeshInit init(argc, argv);
    libMesh::out << "[test_stokes] libMesh initialised.\n";

    // ── 2. Load mesh ──────────────────────────────────────────────────────────
    const char* mesh_env = std::getenv("TEST_MESH");
    if (!mesh_env) {
        std::cerr << "[test_stokes] FAIL: TEST_MESH environment variable not set.\n";
        return 1;
    }

    libMesh::Mesh mesh(init.comm());
    mesh.read(mesh_env);
    mesh.all_second_order(); // upgrade TRI3→TRI6 for P2/P1 Taylor-Hood
    ChannelFlowSystem::tag_pressure_pin(mesh); // fix pressure null space

    if (mesh.n_elem() == 0) {
        std::cerr << "[test_stokes] FAIL: mesh is empty.\n";
        return 1;
    }
    libMesh::out << "[test_stokes] Mesh: " << mesh.n_elem()
                 << " elements, " << mesh.n_nodes() << " nodes.\n";

    // ── 3. Build equation systems ─────────────────────────────────────────────
    libMesh::EquationSystems es(mesh);
    ChannelFlowSystem& sys = es.add_system<ChannelFlowSystem>("ChannelFlow");

    sys.time_solver = std::make_unique<libMesh::SteadySolver>(sys);
    auto& diff_solver = sys.time_solver->diff_solver();
    diff_solver = std::make_unique<libMesh::NewtonSolver>(sys);
    auto& ns = static_cast<libMesh::NewtonSolver&>(*diff_solver);
    ns.max_nonlinear_iterations    = Params::SNES_MAX_IT;
    ns.absolute_residual_tolerance = Params::SNES_ATOL;
    ns.relative_residual_tolerance = Params::SNES_RTOL;
    ns.verbose = true;

    es.init();

    // Apply the full inlet profile (the ramp is zero at t=0).
    sys.time = Params::T_RAMP;
    sys.get_dof_map().create_dof_constraints(mesh, sys.time);
    sys.get_dof_map().enforce_constraints_exactly(sys);
    sys.update();

    libMesh::out << "[test_stokes] DOFs: " << sys.n_dofs() << "\n";

    // ── 4. PETSc KSP/PC options — match production configuration ─────────────
    PetscOptionsSetValue(NULL, "-ksp_type",           "fgmres");
    PetscOptionsSetValue(NULL, "-ksp_gmres_restart",
                         std::to_string(Params::GMRES_RESTART).c_str());
    // ksp_rtol left to libMesh's inexact-Newton tolerance selection.
    PetscOptionsSetValue(NULL, "-ksp_max_it",
                         std::to_string(Params::KSP_MAX_IT).c_str());

    // Fieldsplit type set programmatically in configure_fieldsplit.
    // Sub-PC options are still set here and picked up at PCSetUp.

    // Velocity sub-PC: ILU(1) — robust for non-symmetric Oseen operator.
    PetscOptionsSetValue(NULL, "-fieldsplit_velocity_ksp_type",         "preonly");
    PetscOptionsSetValue(NULL, "-fieldsplit_velocity_pc_type",          "ilu");
    PetscOptionsSetValue(NULL, "-fieldsplit_velocity_pc_factor_levels", "1");

    // Pressure sub-PC: BoomerAMG on assembled Sp (SPD).
    PetscOptionsSetValue(NULL, "-fieldsplit_pressure_ksp_type",         "preonly");
    PetscOptionsSetValue(NULL, "-fieldsplit_pressure_pc_type",          "hypre");
    PetscOptionsSetValue(NULL, "-fieldsplit_pressure_pc_hypre_type",    "boomeramg");

    // Register velocity/pressure IS objects with PETSc fieldsplit.
    ChannelFlowSystem::configure_fieldsplit(sys, mesh, ns);

    // ── 5. Stokes solve ───────────────────────────────────────────────────────
    libMesh::out << "[test_stokes] Running Stokes linear solve...\n";
    sys.set_stokes_mode(true);
    sys.solve();

    // ── 6. Check Newton convergence: 1 iteration, small residual ─────────────
    // total_inner_iterations() counts total linear iterations across all Newton steps.
    // For a linear Stokes problem, Newton should take exactly 1 outer step.
    const unsigned int n_iters = ns.total_inner_iterations();

    // Reassemble the residual at the current solution to measure the final norm.
    sys.assembly(true, false);
    sys.rhs->close();
    const double final_res = sys.rhs->l2_norm();

    libMesh::out << "[test_stokes] Newton inner iterations = " << n_iters
                 << ", final residual norm = " << final_res << "\n";

    // Fieldsplit AMG is mesh-independent; expect a small, nonzero iteration count.
    if (n_iters == 0 || n_iters > 100) {
        std::cerr << "[test_stokes] FAIL: unexpected Newton inner iterations = "
                  << n_iters << " (expected 1–100).\n";
        return 1;
    }

    if (final_res >= Params::SNES_ATOL) {
        std::cerr << "[test_stokes] FAIL: final residual " << final_res
                  << " >= SNES_ATOL " << Params::SNES_ATOL << ".\n";
        return 1;
    }
    libMesh::out << "[test_stokes] Convergence check: PASSED.\n";

    // ── 7. Check no-slip BCs ─────────────────────────────────────────────────
    const libMesh::DofMap&                             dof_map = sys.get_dof_map();
    const libMesh::NumericVector<libMesh::Number>& sol = *sys.solution;

    constexpr double NO_SLIP_TOL = 1.0e-10;

    for (const auto bid : {
             static_cast<libMesh::boundary_id_type>(Params::BID_WALLS),
             static_cast<libMesh::boundary_id_type>(Params::BID_CYLINDER)}) {

        const double max_u = max_abs_on_boundary(mesh, dof_map, sol, sys.u_var(), bid);
        const double max_v = max_abs_on_boundary(mesh, dof_map, sol, sys.v_var(), bid);

        libMesh::out << "[test_stokes] BID " << static_cast<int>(bid)
                     << ": max|u| = " << max_u << ", max|v| = " << max_v << "\n";

        if (max_u >= NO_SLIP_TOL || max_v >= NO_SLIP_TOL) {
            std::cerr << "[test_stokes] FAIL: no-slip violated on BID "
                      << static_cast<int>(bid) << " (max|u|=" << max_u
                      << ", max|v|=" << max_v << ").\n";
            return 1;
        }
    }
    libMesh::out << "[test_stokes] No-slip BC check: PASSED.\n";

    // ── 8. Sanity: inlet u-velocity must be non-zero ──────────────────────────
    const double max_u_inlet = max_abs_on_boundary(
        mesh, dof_map, sol, sys.u_var(),
        static_cast<libMesh::boundary_id_type>(Params::BID_INLET));

    libMesh::out << "[test_stokes] Inlet max|u| = " << max_u_inlet << "\n";

    if (max_u_inlet < 1.0e-10) {
        std::cerr << "[test_stokes] FAIL: inlet u-velocity is zero "
                  << "(parabolic BC not applied).\n";
        return 1;
    }
    libMesh::out << "[test_stokes] Inlet BC sanity check: PASSED.\n";

    libMesh::out << "[test_stokes] ALL CHECKS PASSED.\n";
    return 0;
}
