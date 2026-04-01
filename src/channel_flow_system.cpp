/**
 * @file channel_flow_system.cpp
 * @brief ChannelFlowSystem: DG variable setup, context building, fieldsplit.
 *
 * Weak-form volume assembly is in channel_flow_assembly.cpp.
 * DG face assembly is in dg_face_assembly.cpp.
 */

#include "channel_flow_system.h"
#include "params.h"

#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/dg_fem_context.h"
#include "libmesh/dense_vector.h"
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/sparse_matrix.h"

#include <petscpc.h>
#include <petscis.h>
#include <petscmat.h>

#include <algorithm>
#include <limits>
#include <vector>

// ── ChannelFlowSystem ─────────────────────────────────────────────────────────

ChannelFlowSystem::ChannelFlowSystem(libMesh::EquationSystems& es,
                                     const std::string& name,
                                     unsigned int number)
    : libMesh::FEMSystem(es, name, number),
      _u_var(0), _v_var(0), _p_var(0)
{}

void ChannelFlowSystem::init_data()
{
    // ── Add DG variables (equal-order L2_LAGRANGE) ───────────────────────────
    const auto order = static_cast<libMesh::Order>(Params::DG_ORDER);
    _u_var = this->add_variable("u", order, libMesh::L2_LAGRANGE);
    _v_var = this->add_variable("v", order, libMesh::L2_LAGRANGE);
    _p_var = this->add_variable("p", order, libMesh::L2_LAGRANGE);

    // No Dirichlet BCs — DG enforces BCs weakly through face fluxes.

    // Enable internal-side assembly for DG face integrals.
    this->compute_internal_sides = true;

    // Mark velocity variables as time-evolving (pressure is algebraic).
    this->time_evolving(_u_var, 1);
    this->time_evolving(_v_var, 1);

    this->FEMSystem::init_data();
}

std::unique_ptr<libMesh::DiffContext> ChannelFlowSystem::build_context()
{
    return std::make_unique<libMesh::DGFEMContext>(*this);
}

void ChannelFlowSystem::init_context(libMesh::DiffContext& ctx)
{
    auto& c = libMesh::cast_ref<libMesh::FEMContext&>(ctx);

    // ── Element interior FE data ─────────────────────────────────────────────
    libMesh::FEBase* u_elem_fe = nullptr;
    libMesh::FEBase* p_elem_fe = nullptr;
    c.get_element_fe(_u_var, u_elem_fe);
    c.get_element_fe(_p_var, p_elem_fe);

    u_elem_fe->get_JxW();
    u_elem_fe->get_phi();
    u_elem_fe->get_dphi();
    p_elem_fe->get_phi();
    p_elem_fe->get_dphi();

    // ── Side FE data (current element side) ──────────────────────────────────
    libMesh::FEBase* u_side_fe = nullptr;
    libMesh::FEBase* p_side_fe = nullptr;
    c.get_side_fe(_u_var, u_side_fe);
    c.get_side_fe(_p_var, p_side_fe);

    u_side_fe->get_JxW();
    u_side_fe->get_phi();
    u_side_fe->get_dphi();
    u_side_fe->get_normals();
    u_side_fe->get_xyz();  // needed for BC evaluation at quadrature points
    p_side_fe->get_phi();
    p_side_fe->get_dphi();

    // ── Neighbor side FE data (for DG internal faces) ────────────────────────
    auto& dg_c = libMesh::cast_ref<libMesh::DGFEMContext&>(ctx);
    libMesh::FEBase* u_neigh_fe = nullptr;
    libMesh::FEBase* p_neigh_fe = nullptr;
    dg_c.get_neighbor_side_fe(_u_var, u_neigh_fe);
    dg_c.get_neighbor_side_fe(_p_var, p_neigh_fe);

    u_neigh_fe->get_phi();
    u_neigh_fe->get_dphi();
    p_neigh_fe->get_phi();
    p_neigh_fe->get_dphi();
}

// ── Fieldsplit preconditioner for DG DOFs ─────────────────────────────────────

void ChannelFlowSystem::configure_fieldsplit(ChannelFlowSystem& sys,
                                              libMesh::MeshBase& mesh,
                                              libMesh::NewtonSolver& newton)
{
    newton.reinit();

    auto& ls = libMesh::cast_ref<
        libMesh::PetscLinearSolver<libMesh::Number>&>(newton.get_linear_solver());
    KSP ksp = ls.ksp();

    PC pc;
    KSPGetPC(ksp, &pc);

    // ── Collect velocity and pressure DOF indices (element-local, no sharing) ─
    const libMesh::DofMap& dof_map = sys.get_dof_map();
    std::vector<libMesh::dof_id_type> dofs;
    std::vector<PetscInt> vel_ids, p_ids;

    for (const auto* elem : mesh.active_local_element_ptr_range()) {
        dof_map.dof_indices(elem, dofs, sys.u_var());
        for (auto d : dofs) vel_ids.push_back(static_cast<PetscInt>(d));

        dof_map.dof_indices(elem, dofs, sys.v_var());
        for (auto d : dofs) vel_ids.push_back(static_cast<PetscInt>(d));

        dof_map.dof_indices(elem, dofs, sys.p_var());
        for (auto d : dofs) p_ids.push_back(static_cast<PetscInt>(d));
    }

    // DG DOFs are element-local — no deduplication needed, but sort for PETSc.
    std::sort(vel_ids.begin(), vel_ids.end());
    std::sort(p_ids.begin(), p_ids.end());

    IS vel_is, p_is;
    ISCreateGeneral(PETSC_COMM_WORLD,
                    static_cast<PetscInt>(vel_ids.size()), vel_ids.data(),
                    PETSC_COPY_VALUES, &vel_is);
    ISCreateGeneral(PETSC_COMM_WORLD,
                    static_cast<PetscInt>(p_ids.size()), p_ids.data(),
                    PETSC_COPY_VALUES, &p_is);

    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
    PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_LOWER);
    PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);
    PCFieldSplitSetIS(pc, "velocity", vel_is);
    PCFieldSplitSetIS(pc, "pressure", p_is);

    ISDestroy(&vel_is);
    ISDestroy(&p_is);

    libMesh::out << "[configure_fieldsplit] velocity DOFs: " << vel_ids.size()
                 << "  pressure DOFs: " << p_ids.size() << "\n";
}

void ChannelFlowSystem::setup_pressure_null_space(ChannelFlowSystem& sys,
                                                   libMesh::MeshBase& mesh,
                                                   libMesh::NewtonSolver& /* newton */)
{
    // With DG and no Dirichlet BCs on pressure, the system has a constant
    // pressure null space.  Pin the first local pressure DOF of the first
    // element near the outlet to p=0 via a Dirichlet constraint.
    const libMesh::DofMap& dof_map = sys.get_dof_map();

    // Find the element nearest the outlet-bottom corner
    const libMesh::Elem* pin_elem = nullptr;
    double min_dist2 = std::numeric_limits<double>::max();
    for (const auto* elem : mesh.active_local_element_ptr_range()) {
        const auto centroid = elem->vertex_average();
        const double dx = centroid(0) - Params::CHANNEL_LENGTH;
        const double dy = centroid(1);
        const double d2 = dx * dx + dy * dy;
        if (d2 < min_dist2) {
            min_dist2 = d2;
            pin_elem = elem;
        }
    }

    if (pin_elem) {
        std::vector<libMesh::dof_id_type> p_dofs;
        dof_map.dof_indices(pin_elem, p_dofs, sys.p_var());
        if (!p_dofs.empty()) {
            // Constrain the first pressure DOF of this element to zero
            libMesh::DofConstraintRow empty_row;
            sys.get_dof_map().add_constraint_row(p_dofs[0], empty_row, 0.0, true);
            libMesh::out << "[setup_pressure_null_space] Pinned pressure DOF "
                         << p_dofs[0] << " to zero.\n";
        }
    }
}
