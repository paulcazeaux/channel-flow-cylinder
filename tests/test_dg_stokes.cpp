/**
 * @file test_dg_stokes.cpp
 * @brief DG Phase 4 test: Stokes solve with weakly enforced BCs.
 *
 * Pass criteria:
 *   1. Stokes solve completes (Newton converges in 1 step).
 *   2. Final residual < SNES_ATOL.
 *   3. Velocity on cylinder wall is small (< 0.05; weak BC, not exact zero).
 *   4. Inlet u-velocity is non-zero.
 *   5. Drag force is finite.
 */

#include "channel_flow_system.h"
#include "drag_lift.h"
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
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature_gauss.h"

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

/// @brief Max |var| at quadrature points on boundary sides (DG: element-local interp).
static double max_abs_on_boundary_dg(
    const ChannelFlowSystem& sys,
    const libMesh::MeshBase& mesh,
    unsigned int var,
    libMesh::boundary_id_type bid)
{
    const libMesh::DofMap& dof_map = sys.get_dof_map();
    const libMesh::NumericVector<libMesh::Number>& sol = *sys.current_local_solution;
    const libMesh::FEType fe_type = sys.variable_type(var);

    auto fe = libMesh::FEBase::build(2, fe_type);
    libMesh::QGauss qrule(1, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule(&qrule);
    const auto& phi = fe->get_phi();

    double max_val = 0.0;
    std::vector<libMesh::dof_id_type> dofs;

    for (const auto* elem : mesh.active_local_element_ptr_range()) {
        for (unsigned int s = 0; s < elem->n_sides(); ++s) {
            if (!mesh.get_boundary_info().has_boundary_id(elem, s, bid))
                continue;
            fe->reinit(elem, s);
            dof_map.dof_indices(elem, dofs, var);
            for (unsigned int qp = 0; qp < qrule.n_points(); ++qp) {
                double val = 0.0;
                for (std::size_t j = 0; j < dofs.size(); ++j)
                    val += phi[j][qp] * sol(dofs[j]);
                max_val = std::max(max_val, std::abs(val));
            }
        }
    }
    return max_val;
}

int main(int argc, char** argv)
{
    libMesh::LibMeshInit init(argc, argv);
    libMesh::out << "[test_dg_stokes] libMesh initialised.\n";

    const char* mesh_env = std::getenv("TEST_MESH");
    if (!mesh_env) {
        std::cerr << "[test_dg_stokes] FAIL: TEST_MESH not set.\n";
        return 1;
    }

    libMesh::Mesh mesh(init.comm());
    mesh.read(mesh_env);

    if (mesh.n_elem() == 0) {
        std::cerr << "[test_dg_stokes] FAIL: mesh is empty.\n";
        return 1;
    }
    libMesh::out << "[test_dg_stokes] Mesh: " << mesh.n_elem()
                 << " elements, " << mesh.n_nodes() << " nodes.\n";

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

    // Set time to T_RAMP so inlet profile is fully active.
    sys.time = Params::T_RAMP;

    libMesh::out << "[test_dg_stokes] DOFs: " << sys.n_dofs() << "\n";

    // PETSc options
    PetscOptionsSetValue(NULL, "-ksp_type",           "fgmres");
    PetscOptionsSetValue(NULL, "-ksp_gmres_restart",
                         std::to_string(Params::GMRES_RESTART).c_str());
    PetscOptionsSetValue(NULL, "-ksp_max_it",
                         std::to_string(Params::KSP_MAX_IT).c_str());

    // Use direct solver (LU) for the DG Stokes test to validate assembly.
    PetscOptionsSetValue(NULL, "-pc_type", "lu");

    // Pin pressure DOF to remove null space
    ChannelFlowSystem::setup_pressure_null_space(sys, mesh, ns);

    // Stokes solve
    libMesh::out << "[test_dg_stokes] Running DG Stokes solve...\n";
    sys.set_stokes_mode(true);

    // Debug: assemble and check matrix properties
    sys.assembly(true, true);
    sys.matrix->close();
    sys.rhs->close();
    libMesh::out << "[test_dg_stokes] Matrix assembled. RHS norm = "
                 << sys.rhs->l2_norm() << "\n";

    // Check for zero rows
    {
        const auto& sol = *sys.solution;
        libMesh::out << "[test_dg_stokes] Solution norm (should be 0) = "
                     << sol.l2_norm() << "\n";
    }

    sys.solve();

    // ── Check 1: Newton convergence ──────────────────────────────────────────
    const unsigned int n_iters = ns.total_inner_iterations();
    sys.assembly(true, false);
    sys.rhs->close();
    const double final_res = sys.rhs->l2_norm();

    libMesh::out << "[test_dg_stokes] inner iterations = " << n_iters
                 << ", final residual = " << final_res << "\n";

    if (n_iters == 0 || n_iters > 500) {
        std::cerr << "[test_dg_stokes] FAIL: unexpected iterations = "
                  << n_iters << ".\n";
        return 1;
    }

    if (final_res >= Params::SNES_ATOL) {
        std::cerr << "[test_dg_stokes] FAIL: residual " << final_res
                  << " >= SNES_ATOL " << Params::SNES_ATOL << ".\n";
        return 1;
    }
    libMesh::out << "[test_dg_stokes] Convergence check: PASSED.\n";

    // ── Check 2: weak no-slip on cylinder ────────────────────────────────────
    const auto cyl_bid = static_cast<libMesh::boundary_id_type>(Params::BID_CYLINDER);
    const double max_u_cyl = max_abs_on_boundary_dg(sys, mesh, sys.u_var(), cyl_bid);
    const double max_v_cyl = max_abs_on_boundary_dg(sys, mesh, sys.v_var(), cyl_bid);

    libMesh::out << "[test_dg_stokes] Cylinder: max|u| = " << max_u_cyl
                 << ", max|v| = " << max_v_cyl << "\n";

    // Weak BC tolerance: not exactly zero, but small
    constexpr double WEAK_BC_TOL = 0.05;
    if (max_u_cyl >= WEAK_BC_TOL || max_v_cyl >= WEAK_BC_TOL) {
        std::cerr << "[test_dg_stokes] FAIL: weak no-slip violated.\n";
        return 1;
    }
    libMesh::out << "[test_dg_stokes] Weak no-slip check: PASSED.\n";

    // ── Check 3: inlet velocity non-zero ─────────────────────────────────────
    const auto inlet_bid = static_cast<libMesh::boundary_id_type>(Params::BID_INLET);
    const double max_u_inlet = max_abs_on_boundary_dg(sys, mesh, sys.u_var(), inlet_bid);

    libMesh::out << "[test_dg_stokes] Inlet max|u| = " << max_u_inlet << "\n";

    if (max_u_inlet < 0.1) {
        std::cerr << "[test_dg_stokes] FAIL: inlet u too small.\n";
        return 1;
    }
    libMesh::out << "[test_dg_stokes] Inlet check: PASSED.\n";

    // ── Check 4: drag is finite ──────────────────────────────────────────────
    const auto drag_lift = compute_drag_lift(sys, mesh);
    const double C_D = Params::DRAG_LIFT_NORM * drag_lift.first;
    libMesh::out << "[test_dg_stokes] C_D = " << C_D << "\n";

    if (!std::isfinite(C_D)) {
        std::cerr << "[test_dg_stokes] FAIL: C_D is not finite.\n";
        return 1;
    }
    libMesh::out << "[test_dg_stokes] Drag check: PASSED.\n";

    libMesh::out << "[test_dg_stokes] ALL CHECKS PASSED.\n";
    return 0;
}
