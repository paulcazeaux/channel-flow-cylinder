/**
 * @file dg_assembly.cc
 * @brief DG matrix assembly: A (SIPG viscous + mass), B (pressure coupling),
 *        C (pressure-jump stabilization), M_V, M_p.
 *
 * Uses deal.II's MeshWorker::mesh_loop with cell, interior face, and boundary
 * face workers.  Each matrix is assembled independently for clean block
 * structure.
 */

#include "dg_navier_stokes.h"
#include "dg_params.h"

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_interface_values.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <iostream>
#include <vector>

using namespace dealii;

namespace {

/// @brief Scratch data for DG assembly on one cell + its faces.
struct ScratchData
{
    ScratchData(const Mapping<2>& mapping,
                const FiniteElement<2>& fe,
                const Quadrature<2>& quad,
                const Quadrature<1>& face_quad)
        : fe_values(mapping, fe, quad,
                    update_values | update_gradients |
                    update_JxW_values | update_quadrature_points)
        , fe_face_values(mapping, fe, face_quad,
                         update_values | update_gradients |
                         update_JxW_values | update_normal_vectors |
                         update_quadrature_points)
        , fe_interface_values(mapping, fe, face_quad,
                              update_values | update_gradients |
                              update_JxW_values | update_normal_vectors |
                              update_quadrature_points)
    {}

    ScratchData(const ScratchData& other)
        : fe_values(other.fe_values.get_mapping(),
                    other.fe_values.get_fe(),
                    other.fe_values.get_quadrature(),
                    other.fe_values.get_update_flags())
        , fe_face_values(other.fe_face_values.get_mapping(),
                         other.fe_face_values.get_fe(),
                         other.fe_face_values.get_quadrature(),
                         other.fe_face_values.get_update_flags())
        , fe_interface_values(other.fe_interface_values.get_mapping(),
                              other.fe_interface_values.get_fe(),
                              other.fe_interface_values.get_quadrature(),
                              other.fe_interface_values.get_update_flags())
    {}

    FEValues<2>          fe_values;
    FEFaceValues<2>      fe_face_values;
    FEInterfaceValues<2> fe_interface_values;
};

/// @brief Copy data: local matrices and DOF indices to scatter into global.
struct CopyData
{
    FullMatrix<double>                          cell_matrix;
    std::vector<types::global_dof_index>        local_dof_indices;

    // For interface (face) terms: two sets of DOF indices
    std::vector<FullMatrix<double>>             face_matrices;  // [0]=ee, [1]=en, [2]=ne, [3]=nn
    std::vector<std::vector<types::global_dof_index>> face_dof_indices; // [0]=elem, [1]=neighbor
};

/// @brief Scatter CopyData into a sparse matrix (cell + face contributions).
void copy_to_matrix(SparseMatrix<double>& matrix, const CopyData& data)
{
    // Cell contribution
    if (data.cell_matrix.m() > 0) {
        for (unsigned int i = 0; i < data.local_dof_indices.size(); ++i)
            for (unsigned int j = 0; j < data.local_dof_indices.size(); ++j)
                if (data.cell_matrix(i, j) != 0.0)
                    matrix.add(data.local_dof_indices[i],
                               data.local_dof_indices[j],
                               data.cell_matrix(i, j));
    }

    // Face contributions: one interface matrix per face, each with its own
    // unified DOF index list (element + neighbor DOFs concatenated).
    for (unsigned int b = 0; b < data.face_matrices.size(); ++b) {
        const auto& fm = data.face_matrices[b];
        if (fm.m() == 0)
            continue;
        const auto& dofs = data.face_dof_indices[b];
        for (unsigned int i = 0; i < fm.m(); ++i)
            for (unsigned int j = 0; j < fm.n(); ++j)
                if (fm(i, j) != 0.0)
                    matrix.add(dofs[i], dofs[j], fm(i, j));
    }
}

} // anonymous namespace

// ═══════════════════════════════════════════════════════════════════════════════
//  Assemble A = M_V + gamma*dt*nu*K_DG  (SIPG viscous + mass)
// ═══════════════════════════════════════════════════════════════════════════════

void DGNavierStokes::assemble_implicit_operator()
{
    TimerOutput::Scope t(timer, "assemble_A");

    matrix_A = 0;
    mass_velocity = 0;

    const double nu = Params::NU;
    const double coeff = Params::IMEX_GAMMA * Params::DT * nu;  // γ·dt·ν
    const unsigned int dpc = fe_velocity_scalar.n_dofs_per_cell();

    const QGauss<2> quad(Params::DG_ORDER + 1);
    const QGauss<1> face_quad(Params::DG_ORDER + 1);

    auto cell_worker = [&](const auto& cell, ScratchData& scratch, CopyData& copy) {
        scratch.fe_values.reinit(cell);
        const auto& JxW  = scratch.fe_values.get_JxW_values();
        const unsigned int n_qp = quad.size();

        copy.cell_matrix.reinit(dpc, dpc);
        copy.local_dof_indices.resize(dpc);
        cell->get_dof_indices(copy.local_dof_indices);
        copy.face_matrices.clear();
        copy.face_dof_indices.clear();

        for (unsigned int q = 0; q < n_qp; ++q) {
            for (unsigned int i = 0; i < dpc; ++i)
                for (unsigned int j = 0; j < dpc; ++j) {
                    // Mass: M_V
                    const double mass = scratch.fe_values.shape_value(i, q) *
                                        scratch.fe_values.shape_value(j, q) * JxW[q];
                    // Stiffness: γ·dt·ν · ∫∇φ_i·∇φ_j
                    const double stiff = coeff *
                        scratch.fe_values.shape_grad(i, q) *
                        scratch.fe_values.shape_grad(j, q) * JxW[q];
                    copy.cell_matrix(i, j) += mass + stiff;
                }
        }
    };

    auto face_worker = [&](const auto& cell, unsigned int f,
                           unsigned int sf,
                           const auto& ncell, unsigned int nf,
                           unsigned int nsf,
                           ScratchData& scratch, CopyData& copy) {
        scratch.fe_interface_values.reinit(cell, f, sf, ncell, nf, nsf);
        const auto& iv = scratch.fe_interface_values;
        const auto& JxW = iv.get_JxW_values();
        const auto& normals = iv.get_normal_vectors();
        const unsigned int n_qp = face_quad.size();
        const unsigned int n_dofs = iv.n_current_interface_dofs();

        const double h = cell->face(f)->diameter();
        const double sigma_h = Params::SIGMA_V / h;

        // Resize face matrices: store as single interface matrix
        // then scatter using interface DOF indices
        copy.face_matrices.emplace_back(n_dofs, n_dofs);
        copy.face_dof_indices.push_back(iv.get_interface_dof_indices());

        for (unsigned int q = 0; q < n_qp; ++q) {
            for (unsigned int i = 0; i < n_dofs; ++i)
                for (unsigned int j = 0; j < n_dofs; ++j) {
                    // SIPG: coeff * (-{∇φ·n}[[φ]] - {∇ψ·n}[[ψ]] + σ/h [[φ]][[ψ]])
                    copy.face_matrices.back()(i, j) += coeff * (
                        - iv.average_of_shape_gradients(i, q) * normals[q]
                          * iv.jump_in_shape_values(j, q)
                        - iv.average_of_shape_gradients(j, q) * normals[q]
                          * iv.jump_in_shape_values(i, q)
                        + sigma_h
                          * iv.jump_in_shape_values(i, q)
                          * iv.jump_in_shape_values(j, q)
                    ) * JxW[q];
                }
        }
    };

    auto boundary_worker = [&](const auto& cell, unsigned int f,
                               ScratchData& scratch, CopyData& copy) {
        scratch.fe_face_values.reinit(cell, f);
        const auto& JxW = scratch.fe_face_values.get_JxW_values();
        const auto& normals = scratch.fe_face_values.get_normal_vectors();
        const unsigned int n_qp = face_quad.size();

        const auto bid = cell->face(f)->boundary_id();
        if (bid == Params::BID_OUTLET)
            return;

        const double h = cell->face(f)->diameter();
        const double sigma_h = Params::SIGMA_V / h;

        std::vector<types::global_dof_index> dof_indices(dpc);
        cell->get_dof_indices(dof_indices);

        FullMatrix<double> bdy_matrix(dpc, dpc);
        for (unsigned int q = 0; q < n_qp; ++q) {
            for (unsigned int i = 0; i < dpc; ++i)
                for (unsigned int j = 0; j < dpc; ++j) {
                    const double grad_i_n =
                        scratch.fe_face_values.shape_grad(i, q) * normals[q];
                    const double grad_j_n =
                        scratch.fe_face_values.shape_grad(j, q) * normals[q];
                    const double phi_i = scratch.fe_face_values.shape_value(i, q);
                    const double phi_j = scratch.fe_face_values.shape_value(j, q);

                    bdy_matrix(i, j) += coeff * (
                        - grad_j_n * phi_i     // consistency
                        - grad_i_n * phi_j     // symmetry
                        + sigma_h * phi_i * phi_j  // penalty
                    ) * JxW[q];
                }
        }

        // Add boundary face contribution to global matrix directly
        for (unsigned int i = 0; i < dpc; ++i)
            for (unsigned int j = 0; j < dpc; ++j)
                if (bdy_matrix(i, j) != 0.0)
                    matrix_A.add(dof_indices[i], dof_indices[j], bdy_matrix(i, j));
    };

    ScratchData scratch(mapping, fe_velocity_scalar, quad, face_quad);
    CopyData    copy;

    auto copier = [&](const CopyData& data) {
        copy_to_matrix(matrix_A, data);
    };

    MeshWorker::mesh_loop(dof_handler_velocity.begin_active(),
                          dof_handler_velocity.end(),
                          cell_worker, copier, scratch, copy,
                          MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_own_interior_faces_once |
                          MeshWorker::assemble_boundary_faces,
                          boundary_worker, face_worker);

    // Also assemble M_V separately (needed for Cahouet-Chabard)
    {
        FEValues<2> fe_values(mapping, fe_velocity_scalar, quad,
                              update_values | update_JxW_values);
        FullMatrix<double> local_mass(dpc, dpc);
        std::vector<types::global_dof_index> dof_indices(dpc);

        for (const auto& cell : dof_handler_velocity.active_cell_iterators()) {
            fe_values.reinit(cell);
            local_mass = 0;
            for (unsigned int q = 0; q < quad.size(); ++q)
                for (unsigned int i = 0; i < dpc; ++i)
                    for (unsigned int j = 0; j < dpc; ++j)
                        local_mass(i, j) += fe_values.shape_value(i, q) *
                                            fe_values.shape_value(j, q) *
                                            fe_values.JxW(q);
            cell->get_dof_indices(dof_indices);
            for (unsigned int i = 0; i < dpc; ++i)
                for (unsigned int j = 0; j < dpc; ++j)
                    mass_velocity.add(dof_indices[i], dof_indices[j],
                                      local_mass(i, j));
        }
    }

    std::cout << "Assembled A (nnz=" << matrix_A.n_nonzero_elements()
              << ") and M_V\n";
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Assemble A = nu * K_DG only (steady Stokes, no mass)
// ═══════════════════════════════════════════════════════════════════════════════

void DGNavierStokes::assemble_viscous_only()
{
    TimerOutput::Scope t(timer, "assemble_A_viscous");

    matrix_A = 0;

    const double nu = Params::NU;
    const unsigned int dpc = fe_velocity_scalar.n_dofs_per_cell();

    const QGauss<2> quad(Params::DG_ORDER + 1);
    const QGauss<1> face_quad(Params::DG_ORDER + 1);

    // Cell: ν ∫ ∇φ_i · ∇φ_j
    auto cell_worker = [&](const auto& cell, ScratchData& scratch, CopyData& copy) {
        scratch.fe_values.reinit(cell);
        const auto& JxW = scratch.fe_values.get_JxW_values();

        copy.cell_matrix.reinit(dpc, dpc);
        copy.local_dof_indices.resize(dpc);
        cell->get_dof_indices(copy.local_dof_indices);
        copy.face_matrices.clear();
        copy.face_dof_indices.clear();

        for (unsigned int q = 0; q < quad.size(); ++q)
            for (unsigned int i = 0; i < dpc; ++i)
                for (unsigned int j = 0; j < dpc; ++j)
                    copy.cell_matrix(i, j) += nu *
                        scratch.fe_values.shape_grad(i, q) *
                        scratch.fe_values.shape_grad(j, q) * JxW[q];
    };

    // Interior faces: ν × SIPG (penalty = σ/h, no extra ν)
    auto face_worker = [&](const auto& cell, unsigned int f,
                           unsigned int sf,
                           const auto& ncell, unsigned int nf,
                           unsigned int nsf,
                           ScratchData& scratch, CopyData& copy) {
        scratch.fe_interface_values.reinit(cell, f, sf, ncell, nf, nsf);
        const auto& iv = scratch.fe_interface_values;
        const auto& JxW = iv.get_JxW_values();
        const auto& normals = iv.get_normal_vectors();
        const unsigned int n_dofs = iv.n_current_interface_dofs();

        const double h = cell->face(f)->diameter();
        const double sigma_h = Params::SIGMA_V / h;

        copy.face_matrices.emplace_back(n_dofs, n_dofs);
        copy.face_dof_indices.push_back(iv.get_interface_dof_indices());

        for (unsigned int q = 0; q < face_quad.size(); ++q)
            for (unsigned int i = 0; i < n_dofs; ++i)
                for (unsigned int j = 0; j < n_dofs; ++j)
                    copy.face_matrices.back()(i, j) += nu * (
                        - iv.average_of_shape_gradients(i, q) * normals[q]
                          * iv.jump_in_shape_values(j, q)
                        - iv.average_of_shape_gradients(j, q) * normals[q]
                          * iv.jump_in_shape_values(i, q)
                        + sigma_h
                          * iv.jump_in_shape_values(i, q)
                          * iv.jump_in_shape_values(j, q)
                    ) * JxW[q];
    };

    // Boundary faces: ν × SIPG (non-outlet Dirichlet)
    auto boundary_worker = [&](const auto& cell, unsigned int f,
                               ScratchData& scratch, CopyData&) {
        const auto bid = cell->face(f)->boundary_id();
        if (bid == Params::BID_OUTLET)
            return;

        scratch.fe_face_values.reinit(cell, f);
        const auto& JxW = scratch.fe_face_values.get_JxW_values();

        const double h = cell->face(f)->diameter();
        const double sigma_h = Params::SIGMA_V / h;

        std::vector<types::global_dof_index> dof_indices(dpc);
        cell->get_dof_indices(dof_indices);

        for (unsigned int q = 0; q < face_quad.size(); ++q)
            for (unsigned int i = 0; i < dpc; ++i)
                for (unsigned int j = 0; j < dpc; ++j) {
                    const auto& n = scratch.fe_face_values.normal_vector(q);
                    const double val = nu * (
                        - scratch.fe_face_values.shape_grad(j, q) * n
                          * scratch.fe_face_values.shape_value(i, q)
                        - scratch.fe_face_values.shape_grad(i, q) * n
                          * scratch.fe_face_values.shape_value(j, q)
                        + sigma_h
                          * scratch.fe_face_values.shape_value(i, q)
                          * scratch.fe_face_values.shape_value(j, q)
                    ) * JxW[q];
                    matrix_A.add(dof_indices[i], dof_indices[j], val);
                }
    };

    ScratchData scratch(mapping, fe_velocity_scalar, quad, face_quad);
    CopyData    copy;

    unsigned int total_cells = 0, total_face_blocks = 0;
    auto copier = [&](const CopyData& data) {
        copy_to_matrix(matrix_A, data);
        ++total_cells;
        total_face_blocks += data.face_matrices.size();
    };

    MeshWorker::mesh_loop(dof_handler_velocity.begin_active(),
                          dof_handler_velocity.end(),
                          cell_worker, copier, scratch, copy,
                          MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_own_interior_faces_once |
                          MeshWorker::assemble_boundary_faces,
                          boundary_worker, face_worker);

    std::cout << "Assembled A = ν·K_DG (steady, nnz="
              << matrix_A.n_nonzero_elements()
              << ", cells=" << total_cells
              << ", face_blocks=" << total_face_blocks << ")\n";
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Assemble B (pressure gradient) and B^T (velocity divergence)
// ═══════════════════════════════════════════════════════════════════════════════

void DGNavierStokes::assemble_pressure_coupling()
{
    TimerOutput::Scope t(timer, "assemble_B");

    matrix_B = 0;
    matrix_Bt = 0;

    const unsigned int vel_dpc = fe_velocity_scalar.n_dofs_per_cell();
    const unsigned int pres_dpc = fe_pressure.n_dofs_per_cell();
    const unsigned int n_scalar = dof_handler_velocity.n_dofs();
    const QGauss<2> quad(Params::DG_ORDER + 1);

    // Strong form: B(i,j) = -∫_K ψ_j ∂φ_i/∂x_d (volume only, NO face terms).
    // Inf-sup stability from pressure-jump stabilization C.
    // B^T is the exact transpose.
    {
        FEValues<2> vel_fe(mapping, fe_velocity_scalar, quad,
                           update_gradients | update_JxW_values);
        FEValues<2> pres_fe(mapping, fe_pressure, quad,
                            update_values | update_JxW_values);

        FullMatrix<double> local_Bx(vel_dpc, pres_dpc);
        FullMatrix<double> local_By(vel_dpc, pres_dpc);
        std::vector<types::global_dof_index> vel_dofs(vel_dpc);
        std::vector<types::global_dof_index> pres_dofs(pres_dpc);

        auto vel_cell = dof_handler_velocity.begin_active();
        auto pres_cell = dof_handler_pressure.begin_active();

        for (; vel_cell != dof_handler_velocity.end(); ++vel_cell, ++pres_cell) {
            vel_fe.reinit(vel_cell);
            pres_fe.reinit(pres_cell);

            local_Bx = 0;
            local_By = 0;
            for (unsigned int q = 0; q < quad.size(); ++q)
                for (unsigned int i = 0; i < vel_dpc; ++i)
                    for (unsigned int j = 0; j < pres_dpc; ++j) {
                        const double psi_j = pres_fe.shape_value(j, q);
                        const auto& grad_phi_i = vel_fe.shape_grad(i, q);
                        const double jxw = vel_fe.JxW(q);

                        // -∫ ψ_j ∂φ_i/∂x_d
                        local_Bx(i, j) -= psi_j * grad_phi_i[0] * jxw;
                        local_By(i, j) -= psi_j * grad_phi_i[1] * jxw;
                    }

            vel_cell->get_dof_indices(vel_dofs);
            pres_cell->get_dof_indices(pres_dofs);

            for (unsigned int i = 0; i < vel_dpc; ++i)
                for (unsigned int j = 0; j < pres_dpc; ++j) {
                    if (local_Bx(i, j) != 0.0) {
                        matrix_B.add(vel_dofs[i], pres_dofs[j], local_Bx(i, j));
                        matrix_Bt.add(pres_dofs[j], vel_dofs[i], local_Bx(i, j));
                    }
                    if (local_By(i, j) != 0.0) {
                        matrix_B.add(vel_dofs[i] + n_scalar, pres_dofs[j],
                                     local_By(i, j));
                        matrix_Bt.add(pres_dofs[j], vel_dofs[i] + n_scalar,
                                      local_By(i, j));
                    }
                }
        }
    }

    // No face terms for strong-form B.

    std::cout << "Assembled B (" << matrix_B.m() << "x" << matrix_B.n()
              << ", nnz=" << matrix_B.n_nonzero_elements() << ") and B^T\n";
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Assemble C (pressure-jump stabilization) and M_p (pressure mass)
// ═══════════════════════════════════════════════════════════════════════════════

void DGNavierStokes::assemble_pressure_stabilization()
{
    TimerOutput::Scope t(timer, "assemble_C");

    matrix_C = 0;
    mass_pressure = 0;

    const unsigned int pres_dpc = fe_pressure.n_dofs_per_cell();
    const QGauss<2> quad(Params::DG_ORDER + 1);
    const QGauss<1> face_quad(Params::DG_ORDER + 1);

    // Pressure mass matrix M_p (volume integral)
    {
        FEValues<2> fe_values(mapping, fe_pressure, quad,
                              update_values | update_JxW_values);
        FullMatrix<double> local_mass(pres_dpc, pres_dpc);
        std::vector<types::global_dof_index> dof_indices(pres_dpc);

        for (const auto& cell : dof_handler_pressure.active_cell_iterators()) {
            fe_values.reinit(cell);
            local_mass = 0;
            for (unsigned int q = 0; q < quad.size(); ++q)
                for (unsigned int i = 0; i < pres_dpc; ++i)
                    for (unsigned int j = 0; j < pres_dpc; ++j)
                        local_mass(i, j) += fe_values.shape_value(i, q) *
                                            fe_values.shape_value(j, q) *
                                            fe_values.JxW(q);
            cell->get_dof_indices(dof_indices);
            for (unsigned int i = 0; i < pres_dpc; ++i)
                for (unsigned int j = 0; j < pres_dpc; ++j)
                    mass_pressure.add(dof_indices[i], dof_indices[j],
                                      local_mass(i, j));
        }
    }

    // Pressure-jump stabilization: γ_p · h_F · ∫_F [[p]][[q]] dS
    {
        ScratchData scratch(mapping, fe_pressure, quad, face_quad);
        CopyData    copy;

        auto cell_worker = [](const auto&, ScratchData&, CopyData& c) {
            c.cell_matrix.reinit(0, 0);
            c.face_matrices.clear();
            c.face_dof_indices.clear();
        };

        auto face_worker = [&](const auto& cell, unsigned int f,
                               unsigned int sf,
                               const auto& ncell, unsigned int nf,
                               unsigned int nsf,
                               ScratchData& scratch, CopyData& copy) {
            scratch.fe_interface_values.reinit(cell, f, sf, ncell, nf, nsf);
            const auto& iv = scratch.fe_interface_values;
            const auto& JxW = iv.get_JxW_values();
            const unsigned int n_dofs = iv.n_current_interface_dofs();

            const double h = cell->face(f)->diameter();
            const double gamma_h = Params::GAMMA_P * h;

            copy.face_matrices.emplace_back(n_dofs, n_dofs);
            copy.face_dof_indices.push_back(iv.get_interface_dof_indices());

            for (unsigned int q = 0; q < face_quad.size(); ++q)
                for (unsigned int i = 0; i < n_dofs; ++i)
                    for (unsigned int j = 0; j < n_dofs; ++j)
                        copy.face_matrices.back()(i, j) += gamma_h *
                            iv.jump_in_shape_values(i, q) *
                            iv.jump_in_shape_values(j, q) * JxW[q];
        };

        unsigned int c_cells = 0, c_faces = 0;
        auto copier = [&](const CopyData& data) {
            copy_to_matrix(matrix_C, data);
            ++c_cells;
            c_faces += data.face_matrices.size();
        };

        MeshWorker::mesh_loop(dof_handler_pressure.begin_active(),
                              dof_handler_pressure.end(),
                              cell_worker, copier, scratch, copy,
                              MeshWorker::assemble_own_cells |
                              MeshWorker::assemble_own_interior_faces_once,
                              /*boundary_worker=*/{}, face_worker);

        std::cout << "  C face_blocks=" << c_faces
                  << " (cells=" << c_cells << ")\n";
    }

    std::cout << "Assembled C (nnz=" << matrix_C.n_nonzero_elements()
              << ") and M_p\n";
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Assemble L_p = B^T M_V^{-1} B (discrete pressure Laplacian)
// ═══════════════════════════════════════════════════════════════════════════════

void DGNavierStokes::assemble_pressure_laplacian()
{
    TimerOutput::Scope t(timer, "assemble_Lp");

    matrix_Lp = 0;

    // Assemble a DG SIPG Laplacian on the pressure space as the
    // preconditioner approximation to the Schur complement.
    // This is spectrally equivalent to B^T M_V^{-1} B and has proper
    // inter-element coupling from the face penalty terms.

    const unsigned int pres_dpc = fe_pressure.n_dofs_per_cell();
    const QGauss<2> quad(Params::DG_ORDER + 1);
    const QGauss<1> face_quad(Params::DG_ORDER + 1);

    // Volume: ∫ ∇ψ_i · ∇ψ_j
    {
        FEValues<2> fe_values(mapping, fe_pressure, quad,
                              update_gradients | update_JxW_values);
        FullMatrix<double> local(pres_dpc, pres_dpc);
        std::vector<types::global_dof_index> dof_indices(pres_dpc);

        for (const auto& cell : dof_handler_pressure.active_cell_iterators()) {
            fe_values.reinit(cell);
            local = 0;
            for (unsigned int q = 0; q < quad.size(); ++q)
                for (unsigned int i = 0; i < pres_dpc; ++i)
                    for (unsigned int j = 0; j < pres_dpc; ++j)
                        local(i, j) += fe_values.shape_grad(i, q) *
                                       fe_values.shape_grad(j, q) *
                                       fe_values.JxW(q);
            cell->get_dof_indices(dof_indices);
            for (unsigned int i = 0; i < pres_dpc; ++i)
                for (unsigned int j = 0; j < pres_dpc; ++j)
                    matrix_Lp.add(dof_indices[i], dof_indices[j], local(i, j));
        }
    }

    // Interior faces: SIPG terms
    {
        ScratchData scratch(mapping, fe_pressure, quad, face_quad);
        CopyData    copy;

        auto cell_worker = [](const auto&, ScratchData&, CopyData& c) {
            c.cell_matrix.reinit(0, 0);
            c.face_matrices.clear();
            c.face_dof_indices.clear();
        };

        auto face_worker = [&](const auto& cell, unsigned int f,
                               unsigned int sf,
                               const auto& ncell, unsigned int nf,
                               unsigned int nsf,
                               ScratchData& scratch, CopyData& copy) {
            scratch.fe_interface_values.reinit(cell, f, sf, ncell, nf, nsf);
            const auto& iv = scratch.fe_interface_values;
            const auto& JxW = iv.get_JxW_values();
            const auto& normals = iv.get_normal_vectors();
            const unsigned int n_dofs = iv.n_current_interface_dofs();

            const double h = cell->face(f)->diameter();
            const double penalty = Params::SIGMA_V / h;  // no ν for pressure Laplacian

            copy.face_matrices.emplace_back(n_dofs, n_dofs);
            copy.face_dof_indices.push_back(iv.get_interface_dof_indices());

            for (unsigned int q = 0; q < face_quad.size(); ++q)
                for (unsigned int i = 0; i < n_dofs; ++i)
                    for (unsigned int j = 0; j < n_dofs; ++j)
                        copy.face_matrices.back()(i, j) += (
                            - iv.average_of_shape_gradients(i, q) * normals[q]
                              * iv.jump_in_shape_values(j, q)
                            - iv.average_of_shape_gradients(j, q) * normals[q]
                              * iv.jump_in_shape_values(i, q)
                            + penalty
                              * iv.jump_in_shape_values(i, q)
                              * iv.jump_in_shape_values(j, q)
                        ) * JxW[q];
        };

        // Boundary: Neumann (do nothing) for outlet; Dirichlet penalty for others
        auto boundary_worker = [&](const auto& cell, unsigned int f,
                                   ScratchData& scratch, CopyData&) {
            const auto bid = cell->face(f)->boundary_id();
            if (bid == Params::BID_OUTLET)
                return;

            scratch.fe_face_values.reinit(cell, f);
            const auto& JxW = scratch.fe_face_values.get_JxW_values();
            const auto& normals = scratch.fe_face_values.get_normal_vectors();

            const double h = cell->face(f)->diameter();
            const double penalty = Params::SIGMA_V / h;

            std::vector<types::global_dof_index> dof_indices(pres_dpc);
            cell->get_dof_indices(dof_indices);

            for (unsigned int q = 0; q < face_quad.size(); ++q)
                for (unsigned int i = 0; i < pres_dpc; ++i)
                    for (unsigned int j = 0; j < pres_dpc; ++j) {
                        const double val = (
                            - scratch.fe_face_values.shape_grad(j, q) * normals[q]
                              * scratch.fe_face_values.shape_value(i, q)
                            - scratch.fe_face_values.shape_grad(i, q) * normals[q]
                              * scratch.fe_face_values.shape_value(j, q)
                            + penalty
                              * scratch.fe_face_values.shape_value(i, q)
                              * scratch.fe_face_values.shape_value(j, q)
                        ) * JxW[q];
                        matrix_Lp.add(dof_indices[i], dof_indices[j], val);
                    }
        };

        auto copier = [&](const CopyData& data) {
            copy_to_matrix(matrix_Lp, data);
        };

        MeshWorker::mesh_loop(dof_handler_pressure.begin_active(),
                              dof_handler_pressure.end(),
                              cell_worker, copier, scratch, copy,
                              MeshWorker::assemble_own_cells |
                              MeshWorker::assemble_own_interior_faces_once |
                              MeshWorker::assemble_boundary_faces,
                              boundary_worker, face_worker);
    }

    std::cout << "Assembled L_p SIPG Laplacian (nnz="
              << matrix_Lp.n_nonzero_elements() << ")\n";
}
