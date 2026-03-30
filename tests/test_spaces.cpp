/**
 * @file test_spaces.cpp
 * @brief Phase 2 test: verify FEM space initialisation for Taylor-Hood P2/P1.
 *
 * Pass criteria (from PLAN.md §Phase 2):
 *   1. libMesh initialises without exception.
 *   2. Mesh loads from TEST_MESH (env var); n_elem > 0.
 *   3. P2 velocity DOFs > P1 pressure DOFs.
 *   4. At least one constrained DOF exists on the inlet boundary (BID 1).
 *   5. At least one constrained DOF exists on the cylinder boundary (BID 4).
 *
 * Run via CTest; mesh path passed through the TEST_MESH environment variable.
 */

#include "channel_flow_system.h"
#include "params.h"

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/dof_map.h"
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"

#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>

/** @brief Return true if any DOF on sides with @p bid is constrained. */
static bool boundary_has_constrained_dof(
    const libMesh::MeshBase& mesh,
    const libMesh::DofMap&   dof_map,
    unsigned int             var,
    libMesh::boundary_id_type bid)
{
    for (const auto& elem : mesh.active_local_element_ptr_range()) {
        for (unsigned short s = 0; s < elem->n_sides(); ++s) {
            if (!mesh.get_boundary_info().has_boundary_id(elem, s, bid))
                continue;
            std::vector<libMesh::dof_id_type> dofs;
            dof_map.dof_indices(elem, dofs, var);
            for (libMesh::dof_id_type d : dofs) {
                if (dof_map.is_constrained_dof(d))
                    return true;
            }
        }
    }
    return false;
}

/** @brief Count unique global DOF IDs for a single variable across all elements. */
static std::size_t count_variable_dofs(
    const libMesh::MeshBase& mesh,
    const libMesh::DofMap&   dof_map,
    unsigned int             var)
{
    std::set<libMesh::dof_id_type> ids;
    std::vector<libMesh::dof_id_type> dofs;
    for (const auto& elem : mesh.active_local_element_ptr_range()) {
        dof_map.dof_indices(elem, dofs, var);
        ids.insert(dofs.begin(), dofs.end());
    }
    return ids.size();
}

int main(int argc, char** argv)
{
    // ── 1. libMesh initialisation ─────────────────────────────────────────────
    libMesh::LibMeshInit init(argc, argv);
    libMesh::out << "[test_spaces] libMesh initialised.\n";

    // ── 2. Load mesh ──────────────────────────────────────────────────────────
    const char* mesh_env = std::getenv("TEST_MESH");
    if (!mesh_env) {
        std::cerr << "[test_spaces] FAIL: TEST_MESH environment variable not set.\n";
        return 1;
    }

    libMesh::Mesh mesh(init.comm());
    mesh.read(mesh_env);
    mesh.all_second_order(); // upgrade TRI3→TRI6 for P2/P1 Taylor-Hood

    if (mesh.n_elem() == 0) {
        std::cerr << "[test_spaces] FAIL: mesh is empty.\n";
        return 1;
    }
    libMesh::out << "[test_spaces] Mesh loaded: " << mesh.n_elem()
                 << " elements, " << mesh.n_nodes() << " nodes.\n";

    // ── 3. Build equation systems ─────────────────────────────────────────────
    libMesh::EquationSystems es(mesh);
    ChannelFlowSystem& sys = es.add_system<ChannelFlowSystem>("ChannelFlow");

    // FEMSystem requires a time solver before init.
    sys.time_solver = std::make_unique<libMesh::SteadySolver>(sys);
    sys.time_solver->diff_solver() = std::make_unique<libMesh::NewtonSolver>(sys);

    es.init();
    libMesh::out << "[test_spaces] EquationSystems initialised.\n";

    // ── 4. Check DOF counts ───────────────────────────────────────────────────
    const libMesh::DofMap& dof_map = sys.get_dof_map();

    const std::size_t n_u = count_variable_dofs(mesh, dof_map, sys.u_var());
    const std::size_t n_p = count_variable_dofs(mesh, dof_map, sys.p_var());

    libMesh::out << "[test_spaces] u DOFs = " << n_u
                 << ", p DOFs = " << n_p << "\n";

    if (n_u <= n_p) {
        std::cerr << "[test_spaces] FAIL: expected n_u_dofs > n_p_dofs "
                  << "(got " << n_u << " vs " << n_p << ").\n";
        return 1;
    }

    // ── 5. Check Dirichlet BCs on inlet (BID 1) ───────────────────────────────
    const bool inlet_ok = boundary_has_constrained_dof(
        mesh, dof_map, sys.u_var(),
        static_cast<libMesh::boundary_id_type>(Params::BID_INLET));

    if (!inlet_ok) {
        std::cerr << "[test_spaces] FAIL: no constrained u DOFs on inlet boundary.\n";
        return 1;
    }
    libMesh::out << "[test_spaces] Inlet Dirichlet DOFs: OK.\n";

    // ── 6. Check Dirichlet BCs on cylinder (BID 4) ───────────────────────────
    const bool cyl_ok = boundary_has_constrained_dof(
        mesh, dof_map, sys.u_var(),
        static_cast<libMesh::boundary_id_type>(Params::BID_CYLINDER));

    if (!cyl_ok) {
        std::cerr << "[test_spaces] FAIL: no constrained u DOFs on cylinder boundary.\n";
        return 1;
    }
    libMesh::out << "[test_spaces] Cylinder Dirichlet DOFs: OK.\n";

    libMesh::out << "[test_spaces] ALL CHECKS PASSED.\n";
    return 0;
}
