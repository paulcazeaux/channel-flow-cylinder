/**
 * @file test_dg_spaces.cpp
 * @brief DG Phase 1 test: verify DG space initialisation with L2_LAGRANGE.
 *
 * Pass criteria:
 *   1. libMesh initialises and mesh loads from TEST_MESH.
 *   2. Variables use L2_LAGRANGE family (discontinuous).
 *   3. Total DOFs = n_elem * local_dofs_per_elem * 3 (u, v, p).
 *   4. No constrained DOFs (DG has no inter-element constraints).
 *   5. compute_internal_sides is enabled.
 */

#include "channel_flow_system.h"
#include "params.h"

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_type.h"
#include "libmesh/elem.h"

#include <cstdlib>
#include <iostream>

int main(int argc, char** argv)
{
    libMesh::LibMeshInit init(argc, argv);
    libMesh::out << "[test_dg_spaces] libMesh initialised.\n";

    const char* mesh_env = std::getenv("TEST_MESH");
    if (!mesh_env) {
        std::cerr << "[test_dg_spaces] FAIL: TEST_MESH not set.\n";
        return 1;
    }

    libMesh::Mesh mesh(init.comm());
    mesh.read(mesh_env);
    // No all_second_order() — DG works on TRI3 elements.

    if (mesh.n_elem() == 0) {
        std::cerr << "[test_dg_spaces] FAIL: mesh is empty.\n";
        return 1;
    }
    libMesh::out << "[test_dg_spaces] Mesh: " << mesh.n_elem()
                 << " elements, " << mesh.n_nodes() << " nodes.\n";

    libMesh::EquationSystems es(mesh);
    ChannelFlowSystem& sys = es.add_system<ChannelFlowSystem>("ChannelFlow");

    sys.time_solver = std::make_unique<libMesh::SteadySolver>(sys);
    sys.time_solver->diff_solver() = std::make_unique<libMesh::NewtonSolver>(sys);

    es.init();

    // ── Check 1: FE family is L2_LAGRANGE ────────────────────────────────────
    const auto u_type = sys.variable_type(sys.u_var());
    const auto v_type = sys.variable_type(sys.v_var());
    const auto p_type = sys.variable_type(sys.p_var());

    if (u_type.family != libMesh::L2_LAGRANGE ||
        v_type.family != libMesh::L2_LAGRANGE ||
        p_type.family != libMesh::L2_LAGRANGE) {
        std::cerr << "[test_dg_spaces] FAIL: expected L2_LAGRANGE for all variables.\n";
        return 1;
    }
    libMesh::out << "[test_dg_spaces] FE family check: PASSED (L2_LAGRANGE).\n";

    // ── Check 2: polynomial order matches DG_ORDER ───────────────────────────
    const auto expected_order = static_cast<libMesh::Order>(Params::DG_ORDER);
    if (u_type.order != expected_order || p_type.order != expected_order) {
        std::cerr << "[test_dg_spaces] FAIL: expected order "
                  << Params::DG_ORDER << " for all variables.\n";
        return 1;
    }
    libMesh::out << "[test_dg_spaces] Order check: PASSED (P"
                 << Params::DG_ORDER << ").\n";

    // ── Check 3: DOF count ───────────────────────────────────────────────────
    // For DG on TRI3 with order k: local DOFs per element = (k+1)(k+2)/2
    const unsigned int local_dofs = (Params::DG_ORDER + 1) * (Params::DG_ORDER + 2) / 2;
    const unsigned int expected_total = mesh.n_elem() * local_dofs * 3; // u, v, p
    const unsigned int actual_total = sys.n_dofs();

    libMesh::out << "[test_dg_spaces] DOFs: " << actual_total
                 << " (expected " << expected_total << ")\n";

    if (actual_total != expected_total) {
        std::cerr << "[test_dg_spaces] FAIL: DOF count mismatch.\n";
        return 1;
    }
    libMesh::out << "[test_dg_spaces] DOF count check: PASSED.\n";

    // ── Check 4: no constrained DOFs ─────────────────────────────────────────
    const libMesh::DofMap& dof_map = sys.get_dof_map();
    const unsigned int n_constrained = dof_map.n_constrained_dofs();
    if (n_constrained > 0) {
        std::cerr << "[test_dg_spaces] FAIL: found " << n_constrained
                  << " constrained DOFs (DG should have none).\n";
        return 1;
    }
    libMesh::out << "[test_dg_spaces] No-constraints check: PASSED.\n";

    // ── Check 5: compute_internal_sides is enabled ───────────────────────────
    if (!sys.compute_internal_sides) {
        std::cerr << "[test_dg_spaces] FAIL: compute_internal_sides is false.\n";
        return 1;
    }
    libMesh::out << "[test_dg_spaces] compute_internal_sides check: PASSED.\n";

    libMesh::out << "[test_dg_spaces] ALL CHECKS PASSED.\n";
    return 0;
}
