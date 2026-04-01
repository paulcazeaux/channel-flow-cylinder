/**
 * @file dg_navier_stokes.cc
 * @brief DG Navier-Stokes solver implementation (deal.II).
 *
 * Phase 1: mesh setup, DG space construction, DOF distribution.
 */

#include "dg_navier_stokes.h"
#include "dg_params.h"

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <cmath>
#include <sys/stat.h>

// ── Constructor ──────────────────────────────────────────────────────────────

DGNavierStokes::DGNavierStokes()
    : mapping(Params::DG_ORDER)
    , fe_velocity_scalar(Params::DG_ORDER)
    , fe_pressure(Params::DG_ORDER - 1)  // Q_{k-1} pressure for inf-sup stability
    , dof_handler_velocity(triangulation)
    , dof_handler_pressure(triangulation)
    , n_vel_dofs(0)
    , n_pres_dofs(0)
    , timer(std::cout, dealii::TimerOutput::summary,
            dealii::TimerOutput::wall_times)
{}

// ── Inlet BC ─────────────────────────────────────────────────────────────────

double DGNavierStokes::inlet_u(double y, double t)
{
    const double H = Params::CHANNEL_HEIGHT;
    const double profile = 4.0 * Params::U_MAX * y * (H - y) / (H * H);
    double ramp = 1.0;
    if (t < Params::T_RAMP) {
        const double s = std::sin(M_PI * t / (2.0 * Params::T_RAMP));
        ramp = s * s;
    }
    return profile * ramp;
}

// ── Mesh setup ───────────────────────────────────────────────────────────────

void DGNavierStokes::setup_mesh(unsigned int n_refinements)
{
    dealii::TimerOutput::Scope t(timer, "setup_mesh");

    dealii::GridGenerator::channel_with_cylinder(
        triangulation,
        /*shell_region_width=*/0.03,
        /*n_shells=*/2,
        /*skewness=*/2.0,
        /*colorize=*/true);

    triangulation.refine_global(n_refinements);

    std::cout << "Mesh: " << triangulation.n_active_cells() << " cells, "
              << triangulation.n_vertices() << " vertices"
              << " (" << n_refinements << " global refinements)\n";

    // Count boundary faces per ID
    std::map<dealii::types::boundary_id, unsigned int> bid_counts;
    unsigned int n_boundary_faces = 0;
    for (const auto& cell : triangulation.active_cell_iterators())
        for (unsigned int f = 0; f < cell->n_faces(); ++f)
            if (cell->at_boundary(f)) {
                ++bid_counts[cell->face(f)->boundary_id()];
                ++n_boundary_faces;
            }
    std::cout << "  Total boundary faces: " << n_boundary_faces << "\n";
    for (const auto& [bid, count] : bid_counts)
        std::cout << "  Boundary ID " << bid << ": " << count << " faces\n";

    // CRITICAL: BID 0 is both "inlet" and the default for untagged faces!
    // If refinement creates new boundary faces, they get BID 0 by default,
    // which would be treated as Dirichlet (inlet) instead of their true type.
    // Check if this is happening:
    if (bid_counts.count(0) && bid_counts[0] > 4 * (1u << n_refinements))
        std::cout << "  WARNING: BID 0 has more faces than expected! "
                  << "Possible untagged boundary faces.\n";
}

// ── DOF setup ────────────────────────────────────────────────────────────────

void DGNavierStokes::setup_dofs()
{
    dealii::TimerOutput::Scope t(timer, "setup_dofs");

    // Velocity: one scalar DG DoFHandler, but we track 2 components (u, v)
    // manually.  Each component has n_scalar_dofs DOFs.
    dof_handler_velocity.distribute_dofs(fe_velocity_scalar);
    dof_handler_pressure.distribute_dofs(fe_pressure);

    const unsigned int n_scalar_vel = dof_handler_velocity.n_dofs();
    n_vel_dofs  = 2 * n_scalar_vel;  // u and v components
    n_pres_dofs = dof_handler_pressure.n_dofs();

    std::cout << "DOFs: velocity = " << n_vel_dofs
              << " (" << n_scalar_vel << " per component)"
              << ", pressure = " << n_pres_dofs
              << ", total = " << n_vel_dofs + n_pres_dofs << "\n";

    // Allocate solution vectors
    solution_velocity.reinit(n_vel_dofs);
    solution_pressure.reinit(n_pres_dofs);
    old_velocity.reinit(n_vel_dofs);
}

// ── Sparsity patterns ────────────────────────────────────────────────────────

void DGNavierStokes::setup_sparsity()
{
    dealii::TimerOutput::Scope t(timer, "setup_sparsity");

    const unsigned int n_scalar_vel = dof_handler_velocity.n_dofs();

    // Velocity-velocity sparsity (A, M_V): scalar DG pattern, used for each
    // component block.  With DG, coupling is through face fluxes (element +
    // face neighbors).
    {
        dealii::DynamicSparsityPattern dsp(n_scalar_vel, n_scalar_vel);
        dealii::DoFTools::make_flux_sparsity_pattern(dof_handler_velocity, dsp);
        sparsity_A.copy_from(dsp);
    }

    // Velocity-pressure sparsity (B): size n_vel × n_pres = (2·n_scalar) × n_pres.
    // Rows [0, n_scalar): u-component coupling via ∂/∂x.
    // Rows [n_scalar, 2·n_scalar): v-component coupling via ∂/∂y.
    // Each scalar DOF i couples to pressure DOFs on the same cell and neighbors.
    {
        dealii::DynamicSparsityPattern dsp(n_vel_dofs, n_pres_dofs);
        const unsigned int vel_dpc = fe_velocity_scalar.n_dofs_per_cell();
        const unsigned int pres_dpc = fe_pressure.n_dofs_per_cell();

        std::vector<dealii::types::global_dof_index> vel_indices(vel_dpc);
        std::vector<dealii::types::global_dof_index> pres_indices(pres_dpc);
        std::vector<dealii::types::global_dof_index> neigh_pres_indices(pres_dpc);

        auto vel_cell = dof_handler_velocity.begin_active();
        auto pres_cell = dof_handler_pressure.begin_active();
        for (; vel_cell != dof_handler_velocity.end(); ++vel_cell, ++pres_cell) {
            vel_cell->get_dof_indices(vel_indices);
            pres_cell->get_dof_indices(pres_indices);

            // Both u and v components couple to pressure on same cell
            for (auto vi : vel_indices)
                for (auto pi : pres_indices) {
                    dsp.add(vi, pi);                       // u block
                    dsp.add(vi + n_scalar_vel, pi);        // v block
                }

            // Face coupling: vel on cell ↔ pres on neighbor
            for (unsigned int f = 0; f < vel_cell->n_faces(); ++f) {
                if (vel_cell->at_boundary(f))
                    continue;
                pres_cell->neighbor(f)->get_dof_indices(neigh_pres_indices);
                for (auto vi : vel_indices)
                    for (auto pi : neigh_pres_indices) {
                        dsp.add(vi, pi);                   // u block
                        dsp.add(vi + n_scalar_vel, pi);    // v block
                    }
            }
        }
        sparsity_B.copy_from(dsp);
    }

    // Pressure-pressure sparsity (C, L_p, M_p): DG flux pattern on pressure DOFs.
    {
        dealii::DynamicSparsityPattern dsp(n_pres_dofs, n_pres_dofs);
        dealii::DoFTools::make_flux_sparsity_pattern(dof_handler_pressure, dsp);
        sparsity_C.copy_from(dsp);
        sparsity_Lp.copy_from(dsp);
    }

    // Initialize matrices with sparsity patterns
    matrix_A.reinit(sparsity_A);
    mass_velocity.reinit(sparsity_A);
    matrix_B.reinit(sparsity_B);
    matrix_C.reinit(sparsity_C);
    matrix_Lp.reinit(sparsity_Lp);
    mass_pressure.reinit(sparsity_C);

    // B^T has transposed sparsity of B: size n_pres × n_vel
    {
        dealii::DynamicSparsityPattern dsp_bt(n_pres_dofs, n_vel_dofs);
        for (unsigned int row = 0; row < n_vel_dofs; ++row)
            for (auto it = sparsity_B.begin(row); it != sparsity_B.end(row); ++it)
                dsp_bt.add(it->column(), row);
        sparsity_Bt.copy_from(dsp_bt);
        matrix_Bt.reinit(sparsity_Bt);
    }

    std::cout << "Sparsity: A=" << sparsity_A.n_nonzero_elements()
              << ", B=" << sparsity_B.n_nonzero_elements()
              << ", C=" << sparsity_C.n_nonzero_elements() << "\n";
}

// ── Mass matrix inverses ─────────────────────────────────────────────────────

void DGNavierStokes::setup_mass_inverses()
{
    dealii::TimerOutput::Scope t(timer, "setup_mass_inverses");

    const unsigned int dofs_per_cell = fe_velocity_scalar.n_dofs_per_cell();
    const dealii::QGauss<2> quadrature(Params::DG_ORDER + 1);
    const auto& mapping = this->mapping;

    // Velocity mass inverses
    {
        dealii::FEValues<2> fe_values(mapping, fe_velocity_scalar, quadrature,
                                      dealii::update_values |
                                      dealii::update_JxW_values);

        mv_inv_blocks.resize(triangulation.n_active_cells());
        dealii::FullMatrix<double> local_mass(dofs_per_cell, dofs_per_cell);

        unsigned int cell_idx = 0;
        for (const auto& cell : dof_handler_velocity.active_cell_iterators()) {
            fe_values.reinit(cell);
            local_mass = 0;

            for (unsigned int q = 0; q < quadrature.size(); ++q) {
                const double jxw = fe_values.JxW(q);
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        local_mass(i, j) += fe_values.shape_value(i, q) *
                                            fe_values.shape_value(j, q) * jxw;
            }

            mv_inv_blocks[cell_idx].reinit(dofs_per_cell, dofs_per_cell);
            mv_inv_blocks[cell_idx].invert(local_mass);
            ++cell_idx;
        }
    }

    // Pressure mass inverses
    {
        const unsigned int p_dofs_per_cell = fe_pressure.n_dofs_per_cell();
        dealii::FEValues<2> fe_values(mapping, fe_pressure, quadrature,
                                      dealii::update_values |
                                      dealii::update_JxW_values);

        mp_inv_blocks.resize(triangulation.n_active_cells());
        dealii::FullMatrix<double> local_mass(p_dofs_per_cell, p_dofs_per_cell);

        unsigned int cell_idx = 0;
        for (const auto& cell : dof_handler_pressure.active_cell_iterators()) {
            fe_values.reinit(cell);
            local_mass = 0;

            for (unsigned int q = 0; q < quadrature.size(); ++q) {
                const double jxw = fe_values.JxW(q);
                for (unsigned int i = 0; i < p_dofs_per_cell; ++i)
                    for (unsigned int j = 0; j < p_dofs_per_cell; ++j)
                        local_mass(i, j) += fe_values.shape_value(i, q) *
                                            fe_values.shape_value(j, q) * jxw;
            }

            mp_inv_blocks[cell_idx].reinit(p_dofs_per_cell, p_dofs_per_cell);
            mp_inv_blocks[cell_idx].invert(local_mass);
            ++cell_idx;
        }
    }

    std::cout << "Mass inverses: " << mv_inv_blocks.size()
              << " velocity blocks (" << dofs_per_cell << "x" << dofs_per_cell
              << "), " << mp_inv_blocks.size() << " pressure blocks\n";
}

// ── Output ───────────────────────────────────────────────────────────────────

void DGNavierStokes::output_results(double time, int step) const
{
    const unsigned int n_scalar = dof_handler_velocity.n_dofs();
    dealii::Vector<double> u_component(n_scalar);
    dealii::Vector<double> v_component(n_scalar);
    for (unsigned int i = 0; i < n_scalar; ++i) {
        u_component[i] = solution_velocity[i];
        v_component[i] = solution_velocity[n_scalar + i];
    }

    dealii::DataOut<2> data_out;
    data_out.attach_dof_handler(dof_handler_velocity);
    data_out.add_data_vector(u_component, "u");
    data_out.add_data_vector(v_component, "v");

    // Add pressure on the velocity mesh (interpolated from Q_{k-1} to Q_k)
    dealii::Vector<double> p_on_vel_mesh(n_scalar);
    {
        const unsigned int p_dpc = fe_pressure.n_dofs_per_cell();
        const unsigned int v_dpc = fe_velocity_scalar.n_dofs_per_cell();
        const dealii::QGauss<2> quad(Params::DG_ORDER + 1);
        dealii::FEValues<2> vel_fe(mapping, fe_velocity_scalar, quad,
                                   dealii::update_values);
        dealii::FEValues<2> pres_fe(mapping, fe_pressure, quad,
                                    dealii::update_values | dealii::update_JxW_values);

        std::vector<dealii::types::global_dof_index> v_dofs(v_dpc);
        std::vector<dealii::types::global_dof_index> p_dofs(p_dpc);
        dealii::FullMatrix<double> local_mass(v_dpc, v_dpc);
        dealii::Vector<double> local_rhs(v_dpc), local_sol(v_dpc);

        auto vel_cell = dof_handler_velocity.begin_active();
        auto pres_cell = dof_handler_pressure.begin_active();
        unsigned int cell_idx = 0;
        for (; vel_cell != dof_handler_velocity.end();
             ++vel_cell, ++pres_cell, ++cell_idx) {
            vel_fe.reinit(vel_cell);
            pres_fe.reinit(pres_cell);
            vel_cell->get_dof_indices(v_dofs);
            pres_cell->get_dof_indices(p_dofs);

            // L2 project pressure onto velocity space (element-local)
            local_rhs = 0;
            for (unsigned int q = 0; q < quad.size(); ++q) {
                double p_val = 0;
                for (unsigned int j = 0; j < p_dpc; ++j)
                    p_val += pres_fe.shape_value(j, q) * solution_pressure(p_dofs[j]);
                for (unsigned int i = 0; i < v_dpc; ++i)
                    local_rhs(i) += vel_fe.shape_value(i, q) * p_val * pres_fe.JxW(q);
            }
            mv_inv_blocks[cell_idx].vmult(local_sol, local_rhs);
            for (unsigned int i = 0; i < v_dpc; ++i)
                p_on_vel_mesh(v_dofs[i]) = local_sol(i);
        }
    }
    data_out.add_data_vector(p_on_vel_mesh, "p");

    data_out.build_patches(Params::DG_ORDER);

    const std::string vtu_name = std::string(Params::OUTPUT_FILE)
                                 + "_" + std::to_string(step) + ".vtu";
    std::ofstream vtu_out(vtu_name);
    data_out.write_vtu(vtu_out);

    // Append to PVD collection file
    pvd_entries.push_back({time, vtu_name});
    std::ofstream pvd_out(std::string(Params::OUTPUT_FILE) + ".pvd");
    dealii::DataOutBase::write_pvd_record(pvd_out, pvd_entries);
}

// ── Stub implementations (later phases) ──────────────────────────────────────

void DGNavierStokes::evaluate_advection(
    const dealii::Vector<double>& velocity,
    dealii::Vector<double>& advection_rhs,
    double time)
{
    dealii::TimerOutput::Scope t(timer, "advection");

    // Conservative DG form of (u·∇)u via integration by parts:
    //   Volume:   -∫_K (u ⊗ u) : ∇w dx  (= -∫ u_c u_d ∂w_c/∂x_d)
    //   Face:     +∫_F F̂·n · w ds  (Lax-Friedrichs numerical flux)
    // Output: advection_rhs = assembled ∫(u·∇)u · w in DG sense.
    // The time loop does rhs -= dt * advection_rhs.

    const unsigned int n_scalar = dof_handler_velocity.n_dofs();
    const unsigned int dpc = fe_velocity_scalar.n_dofs_per_cell();
    const dealii::QGauss<2> quad(Params::DG_ORDER + 1);
    const dealii::QGauss<1> face_quad(Params::DG_ORDER + 1);

    advection_rhs.reinit(n_vel_dofs);

    // ── Volume: -∫_K u_c (u · ∇w_c) dx ──────────────────────────────────
    {
        dealii::FEValues<2> fe_values(mapping, fe_velocity_scalar, quad,
                                      dealii::update_values |
                                      dealii::update_gradients |
                                      dealii::update_JxW_values);

        std::vector<dealii::types::global_dof_index> dof_indices(dpc);

        for (const auto& cell : dof_handler_velocity.active_cell_iterators()) {
            fe_values.reinit(cell);
            cell->get_dof_indices(dof_indices);

            for (unsigned int q = 0; q < quad.size(); ++q) {
                double u_val = 0, v_val = 0;
                for (unsigned int j = 0; j < dpc; ++j) {
                    u_val += fe_values.shape_value(j, q) * velocity(dof_indices[j]);
                    v_val += fe_values.shape_value(j, q) * velocity(n_scalar + dof_indices[j]);
                }

                const double jxw = fe_values.JxW(q);

                for (unsigned int i = 0; i < dpc; ++i) {
                    const auto& grad_w = fe_values.shape_grad(i, q);
                    // -∫ u_c (u · ∇w_c) = -(u · ∇w_c) u_c
                    const double u_dot_grad_w = u_val * grad_w[0] + v_val * grad_w[1];
                    advection_rhs(dof_indices[i])             -= u_val * u_dot_grad_w * jxw;
                    advection_rhs(n_scalar + dof_indices[i])  -= v_val * u_dot_grad_w * jxw;
                }
            }
        }
    }

    // ── Faces: +∫_F F̂·n · w ds  (Lax-Friedrichs) ────────────────────────
    {
        dealii::FEFaceValues<2> fe_face(mapping, fe_velocity_scalar, face_quad,
                                        dealii::update_values |
                                        dealii::update_JxW_values |
                                        dealii::update_normal_vectors |
                                        dealii::update_quadrature_points);
        dealii::FEFaceValues<2> fe_face_n(mapping, fe_velocity_scalar, face_quad,
                                          dealii::update_values);

        std::vector<dealii::types::global_dof_index> dofs(dpc), dofs_n(dpc);

        for (const auto& cell : dof_handler_velocity.active_cell_iterators()) {
            cell->get_dof_indices(dofs);

            for (unsigned int f = 0; f < cell->n_faces(); ++f) {
                if (cell->at_boundary(f)) {
                    fe_face.reinit(cell, f);
                    const auto bid = cell->face(f)->boundary_id();

                    for (unsigned int q = 0; q < face_quad.size(); ++q) {
                        const auto& n = fe_face.normal_vector(q);
                        const auto& pt = fe_face.quadrature_point(q);
                        const double jxw = fe_face.JxW(q);

                        double u_int = 0, v_int = 0;
                        for (unsigned int j = 0; j < dpc; ++j) {
                            u_int += fe_face.shape_value(j, q) * velocity(dofs[j]);
                            v_int += fe_face.shape_value(j, q) * velocity(n_scalar + dofs[j]);
                        }

                        double u_ext = 0, v_ext = 0;
                        if (bid == Params::BID_INLET) {
                            u_ext = inlet_u(pt[1], time);
                        } else if (bid == Params::BID_OUTLET) {
                            u_ext = u_int;
                            v_ext = v_int;
                        }

                        const double un_int = u_int * n[0] + v_int * n[1];
                        const double un_ext = u_ext * n[0] + v_ext * n[1];
                        const double lambda = std::max(std::abs(un_int), std::abs(un_ext));

                        // F̂_c · n = ½(u_c^int * un^int + u_c^ext * un^ext)
                        //          + ½λ(u_c^int - u_c^ext)
                        const double flux_u = 0.5 * (u_int * un_int + u_ext * un_ext)
                                            + 0.5 * lambda * (u_int - u_ext);
                        const double flux_v = 0.5 * (v_int * un_int + v_ext * un_ext)
                                            + 0.5 * lambda * (v_int - v_ext);

                        for (unsigned int i = 0; i < dpc; ++i) {
                            const double w_i = fe_face.shape_value(i, q);
                            advection_rhs(dofs[i])            += flux_u * w_i * jxw;
                            advection_rhs(n_scalar + dofs[i]) += flux_v * w_i * jxw;
                        }
                    }
                    continue;
                }

                // Interior face: assemble once per face
                if (!cell->neighbor(f)->is_active() ||
                    cell->active_cell_index() > cell->neighbor(f)->active_cell_index())
                    continue;

                const unsigned int nf = cell->neighbor_of_neighbor(f);
                auto neigh = cell->neighbor(f);

                fe_face.reinit(cell, f);
                fe_face_n.reinit(neigh, nf);
                neigh->get_dof_indices(dofs_n);

                for (unsigned int q = 0; q < face_quad.size(); ++q) {
                    const auto& n = fe_face.normal_vector(q);
                    const double jxw = fe_face.JxW(q);

                    double u_p = 0, v_p = 0;
                    for (unsigned int j = 0; j < dpc; ++j) {
                        u_p += fe_face.shape_value(j, q) * velocity(dofs[j]);
                        v_p += fe_face.shape_value(j, q) * velocity(n_scalar + dofs[j]);
                    }
                    double u_m = 0, v_m = 0;
                    for (unsigned int j = 0; j < dpc; ++j) {
                        u_m += fe_face_n.shape_value(j, q) * velocity(dofs_n[j]);
                        v_m += fe_face_n.shape_value(j, q) * velocity(n_scalar + dofs_n[j]);
                    }

                    const double un_p = u_p * n[0] + v_p * n[1];
                    const double un_m = u_m * n[0] + v_m * n[1];
                    const double lambda = std::max(std::abs(un_p), std::abs(un_m));

                    const double flux_u = 0.5 * (u_p * un_p + u_m * un_m)
                                        + 0.5 * lambda * (u_p - u_m);
                    const double flux_v = 0.5 * (v_p * un_p + v_m * un_m)
                                        + 0.5 * lambda * (v_p - v_m);

                    // + side (cell): flux enters with + sign
                    for (unsigned int i = 0; i < dpc; ++i) {
                        const double w_i = fe_face.shape_value(i, q);
                        advection_rhs(dofs[i])            += flux_u * w_i * jxw;
                        advection_rhs(n_scalar + dofs[i]) += flux_v * w_i * jxw;
                    }
                    // - side (neighbor): flux enters with - sign (outward normal flips)
                    for (unsigned int i = 0; i < dpc; ++i) {
                        const double w_i = fe_face_n.shape_value(i, q);
                        advection_rhs(dofs_n[i])            -= flux_u * w_i * jxw;
                        advection_rhs(n_scalar + dofs_n[i]) -= flux_v * w_i * jxw;
                    }
                }
            }
        }
    }
}

// solve_schur is implemented in dg_schur_solver.cc

std::pair<double, double> DGNavierStokes::compute_drag_lift() const
{
    // TODO: Phase 7 — boundary stress integration
    return {0.0, 0.0};
}

// ── Assemble SIPG BC RHS for velocity ─────────────────────────────────────────

/// @brief Add SIPG weak Dirichlet RHS to velocity vector (u-component only).
/// coeff scales the SIPG terms (= γ·dt·ν for IMEX, or ν for steady).
void DGNavierStokes::assemble_bc_rhs(dealii::Vector<double>& rhs_u,
                                     double coeff, double t) const
{
    const dealii::QGauss<1> face_quad(Params::DG_ORDER + 1);
    const unsigned int dpc = fe_velocity_scalar.n_dofs_per_cell();

    dealii::FEFaceValues<2> fe_face(mapping, fe_velocity_scalar, face_quad,
                                    dealii::update_values |
                                    dealii::update_gradients |
                                    dealii::update_JxW_values |
                                    dealii::update_normal_vectors |
                                    dealii::update_quadrature_points);

    std::vector<dealii::types::global_dof_index> dof_indices(dpc);

    for (const auto& cell : dof_handler_velocity.active_cell_iterators()) {
        for (unsigned int f = 0; f < cell->n_faces(); ++f) {
            if (!cell->at_boundary(f))
                continue;
            const auto bid = cell->face(f)->boundary_id();
            if (bid == Params::BID_OUTLET)
                continue;

            fe_face.reinit(cell, f);
            cell->get_dof_indices(dof_indices);

            const double h = cell->face(f)->diameter();
            const double sigma_h = Params::SIGMA_V / h;

            for (unsigned int q = 0; q < face_quad.size(); ++q) {
                const auto& n = fe_face.normal_vector(q);
                const auto& pt = fe_face.quadrature_point(q);
                const double jxw = fe_face.JxW(q);

                double u_bc = 0.0;
                if (bid == Params::BID_INLET)
                    u_bc = inlet_u(pt[1], t);

                for (unsigned int i = 0; i < dpc; ++i) {
                    const double phi_i = fe_face.shape_value(i, q);
                    const double grad_n = fe_face.shape_grad(i, q) * n;

                    // SIPG BC RHS: coeff * (-symmetry + σ/h penalty)
                    rhs_u(dof_indices[i]) += coeff * (
                        -grad_n * u_bc + sigma_h * u_bc * phi_i
                    ) * jxw;
                }
            }
        }
    }
}

// ── Apply M_V to a velocity vector (element-local) ──────────────────────────

void DGNavierStokes::apply_mass(const dealii::Vector<double>& src,
                                dealii::Vector<double>& dst) const
{
    const unsigned int n_scalar = dof_handler_velocity.n_dofs();
    const unsigned int dpc = fe_velocity_scalar.n_dofs_per_cell();
    std::vector<dealii::types::global_dof_index> dof_indices(dpc);
    dealii::Vector<double> local_src(dpc), local_dst(dpc);

    dst = 0;
    for (unsigned int d = 0; d < 2; ++d) {
        unsigned int cell_idx = 0;
        for (const auto& cell : dof_handler_velocity.active_cell_iterators()) {
            cell->get_dof_indices(dof_indices);
            for (unsigned int i = 0; i < dpc; ++i)
                local_src(i) = src(d * n_scalar + dof_indices[i]);

            // M_V is the inverse of mv_inv — but we stored the inverse.
            // We need M_V * src.  Invert the inverse? No — just use the
            // mass_velocity sparse matrix.
            // Actually, let's use the sparse matrix directly.
            break; // This approach won't work element-locally. Use sparse matvec.
        }
    }
    // Use the sparse mass_velocity matrix instead (one component at a time)
    for (unsigned int d = 0; d < 2; ++d) {
        dealii::Vector<double> src_d(n_scalar), dst_d(n_scalar);
        for (unsigned int i = 0; i < n_scalar; ++i)
            src_d(i) = src(d * n_scalar + i);
        mass_velocity.vmult(dst_d, src_d);
        for (unsigned int i = 0; i < n_scalar; ++i)
            dst(d * n_scalar + i) = dst_d(i);
    }
}

// ── Run ──────────────────────────────────────────────────────────────────────

void DGNavierStokes::run()
{
    std::cout << "DG Navier-Stokes solver (deal.II)\n"
              << "  Re = " << Params::RE
              << ", DG order = " << Params::DG_ORDER
              << ", dt = " << Params::DT << "\n";

    setup_mesh(/*n_refinements=*/0);
    setup_dofs();
    setup_sparsity();
    setup_mass_inverses();

    assemble_implicit_operator();  // A = M_V + γ·dt·ν·K_DG (constant)
    assemble_pressure_coupling();  // B, B^T (strong form, constant)

    // Build element-local L_p^{-1} blocks for Cahouet-Chabard (once)
    {
        dealii::TimerOutput::Scope t(timer, "setup_lp_inv");
        const unsigned int vel_dpc = fe_velocity_scalar.n_dofs_per_cell();
        const unsigned int pres_dpc = fe_pressure.n_dofs_per_cell();
        const dealii::QGauss<2> quad(Params::DG_ORDER + 1);
        dealii::FEValues<2> vel_fe(mapping, fe_velocity_scalar, quad,
                                   dealii::update_gradients | dealii::update_JxW_values);
        dealii::FEValues<2> pres_fe(mapping, fe_pressure, quad,
                                    dealii::update_values | dealii::update_JxW_values);

        lp_inv_blocks.resize(triangulation.n_active_cells());
        dealii::FullMatrix<double> local_B(vel_dpc, pres_dpc);
        dealii::FullMatrix<double> MvinvB(vel_dpc, pres_dpc);
        dealii::FullMatrix<double> local_Lp(pres_dpc, pres_dpc);

        auto vel_cell = dof_handler_velocity.begin_active();
        auto pres_cell = dof_handler_pressure.begin_active();
        unsigned int idx = 0;
        for (; vel_cell != dof_handler_velocity.end();
             ++vel_cell, ++pres_cell, ++idx) {
            vel_fe.reinit(vel_cell);
            pres_fe.reinit(pres_cell);
            local_Lp = 0;
            for (unsigned int dd = 0; dd < 2; ++dd) {
                local_B = 0;
                for (unsigned int q = 0; q < quad.size(); ++q)
                    for (unsigned int i = 0; i < vel_dpc; ++i)
                        for (unsigned int j = 0; j < pres_dpc; ++j)
                            local_B(i, j) -= pres_fe.shape_value(j, q) *
                                vel_fe.shape_grad(i, q)[dd] * vel_fe.JxW(q);
                mv_inv_blocks[idx].mmult(MvinvB, local_B);
                local_B.Tmmult(local_Lp, MvinvB, true);
            }
            lp_inv_blocks[idx].reinit(pres_dpc, pres_dpc);
            lp_inv_blocks[idx].invert(local_Lp);
        }
        std::cout << "L_p^{-1} blocks computed.\n";
    }

    (void)mkdir("results", 0755);
    output_results(0.0, 0);

    const double dt = Params::DT;
    const double coeff = Params::IMEX_GAMMA * dt * Params::NU;
    const unsigned int n_scalar = dof_handler_velocity.n_dofs();

    // ── IMEX-ARK3(2)4L[2]SA (Kennedy & Carpenter 2003) ─────────────────
    // 4 stages, 3rd order, L-stable SDIRK implicit part.
    // Implicit diagonal coefficient γ = IMEX_GAMMA (same for all stages).
    //
    // Explicit Butcher tableau (ã_{ij}):
    //   0    |
    //   c2   | ã21
    //   c3   | ã31  ã32
    //   1    | ã41  ã42  ã43
    //   -----+------------------
    //        | b1   b2   b3   b4
    //
    // Implicit Butcher tableau (a_{ij}):
    //   0    | 0
    //   c2   | a21  γ
    //   c3   | a31  a32  γ
    //   1    | a41  a42  a43  γ
    //   -----+---------------------
    //        | b1   b2   b3   b4
    //
    // The implicit operator A = M_V + γ·dt·ν·K_DG is the SAME for all stages.

    constexpr double g = Params::IMEX_GAMMA;  // 0.4358665215

    // ARK3(2)4L[2]SA coefficients
    constexpr double c2 = 2.0 * g;
    constexpr double c3 = 1.0;  // actually (6g-1)/(4g) but for this tableau c3=1-g ≈ 0.5641

    // Explicit coefficients
    constexpr double ae21 = 2.0 * g;
    constexpr double ae31 = (6.0*g - 1.0) / (4.0*g);
    constexpr double ae32 = (1.0 - 2.0*g) / (4.0*g); // Corrected: should be derived from tableau
    // Actually, let me use the exact ARK3(2)4L[2]SA from the paper:
    // Using the simplified IMEX-SSP3(3,3,2) instead which is simpler:
    // Or just use the standard 3rd-order IMEX from Ascher, Ruuth, Spiteri (1997)

    // Let me use IMEX-Midpoint (2nd order, 1 stage) for simplicity first,
    // then upgrade. Actually, let's just do multi-stage IMEX-Euler (Lie splitting):
    // Each stage: u^{(s)} = u_n + dt * Σ ã_sj N(u^{(j)}) + dt * Σ a_sj L(u^{(j)}, p^{(j)})
    // where L is the implicit Stokes operator.

    // For now: use a simple 2nd-order IMEX (CNAB - Crank-Nicolson / Adams-Bashforth):
    // u_{n+1} = u_n + dt·[-½N(u_n) - ½N(u_{n-1})] + dt·[½L(u_{n+1}) + ½L(u_n)]
    // This requires storing u_{n-1} for the advection extrapolation.
    // But our operator A = M_V + γdt·ν·K_DG has γ=0.5 for Crank-Nicolson.

    // Actually, the simplest robust IMEX: just use IMEX-Euler (1st order) with
    // the correct A that was assembled with IMEX_GAMMA. This is what we already had.
    // Let's keep IMEX-Euler for now and focus on getting vortex shedding.

    int step = 0;
    int output_step = 1;
    for (double t = dt; t <= Params::T_FINAL + 0.5 * dt; t += dt, ++step)
    {
        if (step % 10 == 0)
            std::cout << "Step " << step + 1 << " t=" << t
                      << " |u|=" << solution_velocity.l2_norm()
                      << " max=" << solution_velocity.linfty_norm() << "\n";

        // RHS = M_V * u_n - dt·N(u_n) + γ·dt·ν·(SIPG BC terms)
        dealii::Vector<double> rhs_u(n_vel_dofs);
        dealii::Vector<double> rhs_p(n_pres_dofs);

        apply_mass(solution_velocity, rhs_u);

        dealii::Vector<double> advection(n_vel_dofs);
        evaluate_advection(solution_velocity, advection, t - dt);
        if (!std::isfinite(advection.l2_norm())) {
            std::cout << "  ADVECTION NaN at t=" << t
                      << ", |vel|=" << solution_velocity.l2_norm()
                      << ", max|vel|=" << solution_velocity.linfty_norm() << "\n";
            break;
        }
        rhs_u.add(-dt, advection);

        assemble_bc_rhs(rhs_u, coeff, t);

        if (!std::isfinite(rhs_u.l2_norm())) {
            std::cout << "  RHS NaN at t=" << t << "\n";
            break;
        }

        dealii::Vector<double> sol_u(n_vel_dofs);
        dealii::Vector<double> sol_p(n_pres_dofs);
        solve_schur(rhs_u, rhs_p, sol_u, sol_p);

        solution_velocity = sol_u;
        solution_pressure = sol_p;

        // Check for blowup
        if (!std::isfinite(sol_u.l2_norm())) {
            std::cout << "  BLOWUP at t=" << t << "\n";
            break;
        }

        if ((step + 1) % Params::OUTPUT_INTERVAL == 0 ||
            t >= Params::T_FINAL - 0.5 * dt)
        {
            output_results(t, output_step);
            ++output_step;
            std::cout << "  Output step " << output_step
                      << ", |u|=" << sol_u.l2_norm()
                      << ", |p|=" << sol_p.l2_norm()
                      << ", max|u|=" << sol_u.linfty_norm() << "\n";
        }
    }

    std::cout << "\nSimulation complete. " << output_step
              << " snapshots written.\n";
}
