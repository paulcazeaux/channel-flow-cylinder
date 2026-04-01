/**
 * @file dg_face_assembly.cpp
 * @brief DG face integrals: SIPG viscous, upwind advection, pressure coupling.
 *
 * Implements side_time_derivative() and side_constraint() for both interior
 * faces (DG inter-element fluxes) and boundary faces (weak Dirichlet BCs).
 *
 * Neighbor solution values are extracted from the global solution vector
 * using DGFEMContext::get_neighbor_dof_indices().
 */

#include "channel_flow_system.h"
#include "params.h"

#include "libmesh/dg_fem_context.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/point.h"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace libMesh;

namespace {

/// @brief Parabolic inlet u-velocity with smooth ramp.
double inlet_u(double y, double t)
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

/// @brief Compute face diameter h_F for a side of an element.
double face_diameter(const Elem& elem, unsigned int side)
{
    auto side_ptr = elem.build_side_ptr(side);
    return side_ptr->hmax();
}

/// @brief Extract neighbor solution coefficients for a variable from global vector.
void get_neighbor_solution(const DGFEMContext& c,
                           const NumericVector<Number>& sol,
                           unsigned int var,
                           std::vector<Number>& vals)
{
    const auto& dof_ids = c.get_neighbor_dof_indices(var);
    vals.resize(dof_ids.size());
    for (std::size_t j = 0; j < dof_ids.size(); ++j)
        vals[j] = sol(dof_ids[j]);
}

} // anonymous namespace

// ═══════════════════════════════════════════════════════════════════════════════
//  side_time_derivative: momentum face terms
// ═══════════════════════════════════════════════════════════════════════════════

bool ChannelFlowSystem::side_time_derivative(bool request_jacobian,
                                              DiffContext& ctx)
{
    auto& c = cast_ref<DGFEMContext&>(ctx);
    const Elem& elem = c.get_elem();
    const unsigned int side = c.get_side();

    // FEMSystem does not automatically set up the DG neighbor context.
    // We must detect interior sides and call neighbor_side_fe_reinit() manually.
    const Elem* neighbor = elem.neighbor_ptr(side);
    if (neighbor) {
        // Interior side — only assemble once per face (when elem id < neighbor id)
        if (elem.id() >= neighbor->id())
            return request_jacobian;
        c.set_neighbor(*neighbor);
        c.neighbor_side_fe_reinit();
    }
    // After this point, c.dg_terms_are_active() is true for interior sides

    // ── Side FE data (current element) ───────────────────────────────────────
    FEBase* u_fe = nullptr;
    FEBase* p_fe = nullptr;
    c.get_side_fe(_u_var, u_fe);
    c.get_side_fe(_p_var, p_fe);

    const auto& JxW    = u_fe->get_JxW();
    const auto& phi    = u_fe->get_phi();
    const auto& dphi   = u_fe->get_dphi();
    const auto& normal = u_fe->get_normals();
    const auto& psi    = p_fe->get_phi();

    const std::size_t n_u = c.get_dof_indices(_u_var).size();
    const std::size_t n_p = c.get_dof_indices(_p_var).size();
    const unsigned int n_qp = c.get_side_qrule().n_points();

    const double nu = Params::NU;
    const double h_F = face_diameter(elem, side);
    const double penalty = Params::SIGMA_V * nu / h_F;
    const Real sol_deriv = c.get_elem_solution_derivative();

    if (c.dg_terms_are_active()) {
        // ─── Interior face: DG inter-element fluxes ──────────────────────────

        FEBase* u_neigh_fe = nullptr;
        FEBase* p_neigh_fe = nullptr;
        c.get_neighbor_side_fe(_u_var, u_neigh_fe);
        c.get_neighbor_side_fe(_p_var, p_neigh_fe);

        const auto& phi_n  = u_neigh_fe->get_phi();
        const auto& dphi_n = u_neigh_fe->get_dphi();
        const auto& psi_n  = p_neigh_fe->get_phi();

        const std::size_t n_u_n = c.get_neighbor_dof_indices(_u_var).size();
        const std::size_t n_p_n = c.get_neighbor_dof_indices(_p_var).size();

        // Extract neighbor solution from global vector
        const auto& sol = *this->current_local_solution;
        std::vector<Number> u_n_sol, v_n_sol, p_n_sol;
        get_neighbor_solution(c, sol, _u_var, u_n_sol);
        get_neighbor_solution(c, sol, _v_var, v_n_sol);
        get_neighbor_solution(c, sol, _p_var, p_n_sol);

        auto& Fu  = c.get_elem_residual(_u_var);
        auto& Fv  = c.get_elem_residual(_v_var);
        auto& Fu_n = c.get_neighbor_residual(_u_var);
        auto& Fv_n = c.get_neighbor_residual(_v_var);

        DenseSubMatrix<Number>* Kuu_ee = nullptr, *Kvv_ee = nullptr;
        DenseSubMatrix<Number>* Kuu_en = nullptr, *Kvv_en = nullptr;
        DenseSubMatrix<Number>* Kuu_ne = nullptr, *Kvv_ne = nullptr;
        DenseSubMatrix<Number>* Kuu_nn = nullptr, *Kvv_nn = nullptr;
        DenseSubMatrix<Number>* Kup_ee = nullptr, *Kvp_ee = nullptr;
        DenseSubMatrix<Number>* Kup_en = nullptr, *Kvp_en = nullptr;
        DenseSubMatrix<Number>* Kup_ne = nullptr, *Kvp_ne = nullptr;
        DenseSubMatrix<Number>* Kup_nn = nullptr, *Kvp_nn = nullptr;
        if (request_jacobian) {
            Kuu_ee = &c.get_elem_elem_jacobian(_u_var, _u_var);
            Kvv_ee = &c.get_elem_elem_jacobian(_v_var, _v_var);
            Kuu_en = &c.get_elem_neighbor_jacobian(_u_var, _u_var);
            Kvv_en = &c.get_elem_neighbor_jacobian(_v_var, _v_var);
            Kuu_ne = &c.get_neighbor_elem_jacobian(_u_var, _u_var);
            Kvv_ne = &c.get_neighbor_elem_jacobian(_v_var, _v_var);
            Kuu_nn = &c.get_neighbor_neighbor_jacobian(_u_var, _u_var);
            Kvv_nn = &c.get_neighbor_neighbor_jacobian(_v_var, _v_var);
            Kup_ee = &c.get_elem_elem_jacobian(_u_var, _p_var);
            Kvp_ee = &c.get_elem_elem_jacobian(_v_var, _p_var);
            Kup_en = &c.get_elem_neighbor_jacobian(_u_var, _p_var);
            Kvp_en = &c.get_elem_neighbor_jacobian(_v_var, _p_var);
            Kup_ne = &c.get_neighbor_elem_jacobian(_u_var, _p_var);
            Kvp_ne = &c.get_neighbor_elem_jacobian(_v_var, _p_var);
            Kup_nn = &c.get_neighbor_neighbor_jacobian(_u_var, _p_var);
            Kvp_nn = &c.get_neighbor_neighbor_jacobian(_v_var, _p_var);
        }

        for (unsigned int qp = 0; qp < n_qp; ++qp) {
            const Real jxw = JxW[qp];
            const Point& n = normal[qp];

            // Interpolate element-side solution
            Number u_e = 0, v_e = 0, p_e = 0;
            Gradient grad_u_e, grad_v_e;
            for (std::size_t j = 0; j < n_u; ++j) {
                u_e += phi[j][qp] * c.get_elem_solution(_u_var)(j);
                v_e += phi[j][qp] * c.get_elem_solution(_v_var)(j);
                grad_u_e += dphi[j][qp] * c.get_elem_solution(_u_var)(j);
                grad_v_e += dphi[j][qp] * c.get_elem_solution(_v_var)(j);
            }
            for (std::size_t j = 0; j < n_p; ++j)
                p_e += psi[j][qp] * c.get_elem_solution(_p_var)(j);

            // Interpolate neighbor-side solution
            Number u_nb = 0, v_nb = 0, p_nb = 0;
            Gradient grad_u_nb, grad_v_nb;
            for (std::size_t j = 0; j < n_u_n; ++j) {
                u_nb += phi_n[j][qp] * u_n_sol[j];
                v_nb += phi_n[j][qp] * v_n_sol[j];
                grad_u_nb += dphi_n[j][qp] * u_n_sol[j];
                grad_v_nb += dphi_n[j][qp] * v_n_sol[j];
            }
            for (std::size_t j = 0; j < n_p_n; ++j)
                p_nb += psi_n[j][qp] * p_n_sol[j];

            // Jumps and averages
            const Number jump_u = u_e - u_nb;
            const Number jump_v = v_e - v_nb;
            const Gradient avg_grad_u = 0.5 * (grad_u_e + grad_u_nb);
            const Gradient avg_grad_v = 0.5 * (grad_v_e + grad_v_nb);
            const Number avg_p = 0.5 * (p_e + p_nb);

            // SIPG: -avg(grad_u).n + penalty * jump_u
            const Number sipg_u = -nu * (avg_grad_u * n) + penalty * jump_u;
            const Number sipg_v = -nu * (avg_grad_v * n) + penalty * jump_v;

            // Lax-Friedrichs advective flux
            const Number u_avg = 0.5 * (u_e + u_nb);
            const Number v_avg = 0.5 * (v_e + v_nb);
            const Number un = u_avg * n(0) + v_avg * n(1);
            const Number lambda = std::abs(un);
            Number adv_u = 0, adv_v = 0;
            if (!_stokes_mode) {
                adv_u = un * u_avg + 0.5 * lambda * jump_u;
                adv_v = un * v_avg + 0.5 * lambda * jump_v;
            }

            // Assemble elem residual
            for (std::size_t i = 0; i < n_u; ++i) {
                const Real phi_i = phi[i][qp];
                const RealGradient& dphi_i = dphi[i][qp];

                Fu(i) += jxw * (sipg_u * phi_i
                    - nu * (dphi_i * n) * jump_u * 0.5  // symmetry
                    + avg_p * n(0) * phi_i);             // pressure
                Fv(i) += jxw * (sipg_v * phi_i
                    - nu * (dphi_i * n) * jump_v * 0.5
                    + avg_p * n(1) * phi_i);

                if (!_stokes_mode) {
                    Fu(i) += jxw * adv_u * phi_i;
                    Fv(i) += jxw * adv_v * phi_i;
                }
            }
            // Assemble neighbor residual
            for (std::size_t i = 0; i < n_u_n; ++i) {
                const Real phi_ni = phi_n[i][qp];
                const RealGradient& dphi_ni = dphi_n[i][qp];

                Fu_n(i) -= jxw * (sipg_u * phi_ni
                    + nu * (dphi_ni * n) * jump_u * 0.5  // symmetry (note sign)
                    - avg_p * n(0) * phi_ni);             // pressure
                Fv_n(i) -= jxw * (sipg_v * phi_ni
                    + nu * (dphi_ni * n) * jump_v * 0.5
                    - avg_p * n(1) * phi_ni);

                if (!_stokes_mode) {
                    Fu_n(i) -= jxw * adv_u * phi_ni;
                    Fv_n(i) -= jxw * adv_v * phi_ni;
                }
            }

            // Jacobian (SIPG viscous + pressure; advection Jacobian omitted for now)
            if (request_jacobian) {
                for (std::size_t i = 0; i < n_u; ++i) {
                    const Real phi_i = phi[i][qp];
                    const Real sym_i = nu * (dphi[i][qp] * n) * 0.5;
                    for (std::size_t j = 0; j < n_u; ++j) {
                        const Real v1 = jxw * sol_deriv * (
                            -nu * 0.5 * (dphi[j][qp] * n) * phi_i
                            + penalty * phi[j][qp] * phi_i
                            - sym_i * phi[j][qp]);
                        (*Kuu_ee)(i,j) += v1;
                        (*Kvv_ee)(i,j) += v1;
                    }
                    for (std::size_t j = 0; j < n_u_n; ++j) {
                        const Real v1 = jxw * sol_deriv * (
                            nu * 0.5 * (dphi_n[j][qp] * n) * phi_i
                            - penalty * phi_n[j][qp] * phi_i
                            + sym_i * phi_n[j][qp]);
                        (*Kuu_en)(i,j) += v1;
                        (*Kvv_en)(i,j) += v1;
                    }
                    for (std::size_t j = 0; j < n_p; ++j) {
                        (*Kup_ee)(i,j) += jxw * sol_deriv * 0.5 * psi[j][qp] * n(0) * phi_i;
                        (*Kvp_ee)(i,j) += jxw * sol_deriv * 0.5 * psi[j][qp] * n(1) * phi_i;
                    }
                    for (std::size_t j = 0; j < n_p_n; ++j) {
                        (*Kup_en)(i,j) += jxw * sol_deriv * 0.5 * psi_n[j][qp] * n(0) * phi_i;
                        (*Kvp_en)(i,j) += jxw * sol_deriv * 0.5 * psi_n[j][qp] * n(1) * phi_i;
                    }
                }
                for (std::size_t i = 0; i < n_u_n; ++i) {
                    const Real phi_ni = phi_n[i][qp];
                    const Real sym_ni = nu * (dphi_n[i][qp] * n) * 0.5;
                    for (std::size_t j = 0; j < n_u; ++j) {
                        const Real v1 = jxw * sol_deriv * (
                            nu * 0.5 * (dphi[j][qp] * n) * phi_ni
                            - penalty * phi[j][qp] * phi_ni
                            - sym_ni * phi[j][qp]);
                        (*Kuu_ne)(i,j) += v1;
                        (*Kvv_ne)(i,j) += v1;
                    }
                    for (std::size_t j = 0; j < n_u_n; ++j) {
                        const Real v1 = jxw * sol_deriv * (
                            -nu * 0.5 * (dphi_n[j][qp] * n) * phi_ni
                            + penalty * phi_n[j][qp] * phi_ni
                            + sym_ni * phi_n[j][qp]);
                        (*Kuu_nn)(i,j) += v1;
                        (*Kvv_nn)(i,j) += v1;
                    }
                    for (std::size_t j = 0; j < n_p; ++j) {
                        (*Kup_ne)(i,j) -= jxw * sol_deriv * 0.5 * psi[j][qp] * n(0) * phi_ni;
                        (*Kvp_ne)(i,j) -= jxw * sol_deriv * 0.5 * psi[j][qp] * n(1) * phi_ni;
                    }
                    for (std::size_t j = 0; j < n_p_n; ++j) {
                        (*Kup_nn)(i,j) -= jxw * sol_deriv * 0.5 * psi_n[j][qp] * n(0) * phi_ni;
                        (*Kvp_nn)(i,j) -= jxw * sol_deriv * 0.5 * psi_n[j][qp] * n(1) * phi_ni;
                    }
                }
            }
        }
    } else {
        // ─── Boundary face: weakly enforced BCs ─────────────────────────────
        auto& Fu = c.get_elem_residual(_u_var);
        auto& Fv = c.get_elem_residual(_v_var);

        DenseSubMatrix<Number>* Kuu = nullptr, *Kvv = nullptr;
        DenseSubMatrix<Number>* Kup = nullptr, *Kvp = nullptr;
        if (request_jacobian) {
            Kuu = &c.get_elem_jacobian(_u_var, _u_var);
            Kvv = &c.get_elem_jacobian(_v_var, _v_var);
            Kup = &c.get_elem_jacobian(_u_var, _p_var);
            Kvp = &c.get_elem_jacobian(_v_var, _p_var);
        }

        const BoundaryInfo& binfo = this->get_mesh().get_boundary_info();
        const auto bid_s = static_cast<unsigned short>(side);

        // Determine which boundary this side belongs to
        int bid_int = -1;
        for (int b : {Params::BID_INLET, Params::BID_OUTLET,
                      Params::BID_WALLS, Params::BID_CYLINDER})
            if (binfo.has_boundary_id(&elem, bid_s, static_cast<boundary_id_type>(b)))
                { bid_int = b; break; }
        if (bid_int < 0) return request_jacobian;

        // Outlet: stress-free Neumann → do nothing
        if (bid_int == Params::BID_OUTLET) return request_jacobian;

        for (unsigned int qp = 0; qp < n_qp; ++qp) {
            const Real jxw = JxW[qp];
            const Point& n = normal[qp];
            const Point& xyz = u_fe->get_xyz()[qp];

            Number u_h = 0, v_h = 0, p_h = 0;
            Gradient grad_u_h, grad_v_h;
            for (std::size_t j = 0; j < n_u; ++j) {
                u_h += phi[j][qp] * c.get_elem_solution(_u_var)(j);
                v_h += phi[j][qp] * c.get_elem_solution(_v_var)(j);
                grad_u_h += dphi[j][qp] * c.get_elem_solution(_u_var)(j);
                grad_v_h += dphi[j][qp] * c.get_elem_solution(_v_var)(j);
            }
            for (std::size_t j = 0; j < n_p; ++j)
                p_h += psi[j][qp] * c.get_elem_solution(_p_var)(j);

            double u_bc = 0.0, v_bc = 0.0;
            if (bid_int == Params::BID_INLET)
                u_bc = inlet_u(xyz(1), this->time);

            const Number jump_u = u_h - u_bc;
            const Number jump_v = v_h - v_bc;

            for (std::size_t i = 0; i < n_u; ++i) {
                const Real phi_i = phi[i][qp];
                const RealGradient& dphi_i = dphi[i][qp];

                // SIPG: consistency + symmetry + penalty
                Fu(i) -= jxw * nu * (grad_u_h * n) * phi_i;
                Fu(i) -= jxw * nu * (dphi_i * n) * jump_u;
                Fu(i) += jxw * penalty * jump_u * phi_i;
                Fv(i) -= jxw * nu * (grad_v_h * n) * phi_i;
                Fv(i) -= jxw * nu * (dphi_i * n) * jump_v;
                Fv(i) += jxw * penalty * jump_v * phi_i;

                // Pressure: p * (w . n)
                Fu(i) += jxw * p_h * n(0) * phi_i;
                Fv(i) += jxw * p_h * n(1) * phi_i;

                // Advection (upwind)
                if (!_stokes_mode) {
                    const Number un = u_h * n(0) + v_h * n(1);
                    if (un >= 0) {
                        Fu(i) += jxw * un * u_h * phi_i;
                        Fv(i) += jxw * un * v_h * phi_i;
                    } else {
                        Fu(i) += jxw * un * u_bc * phi_i;
                        Fv(i) += jxw * un * v_bc * phi_i;
                    }
                }

                if (request_jacobian) {
                    for (std::size_t j = 0; j < n_u; ++j) {
                        const Real val = jxw * sol_deriv * (
                            -nu * (dphi[j][qp] * n) * phi_i
                            - nu * (dphi_i * n) * phi[j][qp]
                            + penalty * phi[j][qp] * phi_i);
                        (*Kuu)(i,j) += val;
                        (*Kvv)(i,j) += val;
                    }
                    for (std::size_t j = 0; j < n_p; ++j) {
                        (*Kup)(i,j) += jxw * sol_deriv * psi[j][qp] * n(0) * phi_i;
                        (*Kvp)(i,j) += jxw * sol_deriv * psi[j][qp] * n(1) * phi_i;
                    }
                }
            }
        }
    }
    return request_jacobian;
}

// ═══════════════════════════════════════════════════════════════════════════════
//  side_constraint: continuity face terms + pressure-jump stabilisation
// ═══════════════════════════════════════════════════════════════════════════════

bool ChannelFlowSystem::side_constraint(bool request_jacobian,
                                         DiffContext& ctx)
{
    auto& c = cast_ref<DGFEMContext&>(ctx);
    const Elem& elem = c.get_elem();
    const unsigned int side = c.get_side();

    // Manual DG neighbor setup (same as in side_time_derivative)
    const Elem* neighbor = elem.neighbor_ptr(side);
    if (neighbor) {
        if (elem.id() >= neighbor->id())
            return request_jacobian;
        c.set_neighbor(*neighbor);
        c.neighbor_side_fe_reinit();
    }

    FEBase* u_fe = nullptr;
    FEBase* p_fe = nullptr;
    c.get_side_fe(_u_var, u_fe);
    c.get_side_fe(_p_var, p_fe);

    const auto& JxW    = u_fe->get_JxW();
    const auto& phi    = u_fe->get_phi();
    const auto& normal = u_fe->get_normals();
    const auto& psi    = p_fe->get_phi();

    const std::size_t n_u = c.get_dof_indices(_u_var).size();
    const std::size_t n_p = c.get_dof_indices(_p_var).size();
    const unsigned int n_qp = c.get_side_qrule().n_points();
    const double h_F = face_diameter(elem, side);

    if (c.dg_terms_are_active()) {
        // ─── Interior face ───────────────────────────────────────────────────

        FEBase* u_neigh_fe = nullptr;
        FEBase* p_neigh_fe = nullptr;
        c.get_neighbor_side_fe(_u_var, u_neigh_fe);
        c.get_neighbor_side_fe(_p_var, p_neigh_fe);

        const auto& phi_n = u_neigh_fe->get_phi();
        const auto& psi_n = p_neigh_fe->get_phi();

        const std::size_t n_u_n = c.get_neighbor_dof_indices(_u_var).size();
        const std::size_t n_p_n = c.get_neighbor_dof_indices(_p_var).size();

        const auto& sol = *this->current_local_solution;
        std::vector<Number> u_n_sol, v_n_sol, p_n_sol;
        get_neighbor_solution(c, sol, _u_var, u_n_sol);
        get_neighbor_solution(c, sol, _v_var, v_n_sol);
        get_neighbor_solution(c, sol, _p_var, p_n_sol);

        auto& Fp   = c.get_elem_residual(_p_var);
        auto& Fp_n = c.get_neighbor_residual(_p_var);

        DenseSubMatrix<Number>* Kpu_ee = nullptr, *Kpv_ee = nullptr;
        DenseSubMatrix<Number>* Kpu_en = nullptr, *Kpv_en = nullptr;
        DenseSubMatrix<Number>* Kpu_ne = nullptr, *Kpv_ne = nullptr;
        DenseSubMatrix<Number>* Kpu_nn = nullptr, *Kpv_nn = nullptr;
        DenseSubMatrix<Number>* Kpp_ee = nullptr, *Kpp_en = nullptr;
        DenseSubMatrix<Number>* Kpp_ne = nullptr, *Kpp_nn = nullptr;
        if (request_jacobian) {
            Kpu_ee = &c.get_elem_elem_jacobian(_p_var, _u_var);
            Kpv_ee = &c.get_elem_elem_jacobian(_p_var, _v_var);
            Kpu_en = &c.get_elem_neighbor_jacobian(_p_var, _u_var);
            Kpv_en = &c.get_elem_neighbor_jacobian(_p_var, _v_var);
            Kpu_ne = &c.get_neighbor_elem_jacobian(_p_var, _u_var);
            Kpv_ne = &c.get_neighbor_elem_jacobian(_p_var, _v_var);
            Kpu_nn = &c.get_neighbor_neighbor_jacobian(_p_var, _u_var);
            Kpv_nn = &c.get_neighbor_neighbor_jacobian(_p_var, _v_var);
            Kpp_ee = &c.get_elem_elem_jacobian(_p_var, _p_var);
            Kpp_en = &c.get_elem_neighbor_jacobian(_p_var, _p_var);
            Kpp_ne = &c.get_neighbor_elem_jacobian(_p_var, _p_var);
            Kpp_nn = &c.get_neighbor_neighbor_jacobian(_p_var, _p_var);
        }

        for (unsigned int qp = 0; qp < n_qp; ++qp) {
            const Real jxw = JxW[qp];
            const Point& n = normal[qp];

            // Interpolate velocities
            Number u_e = 0, v_e = 0;
            for (std::size_t j = 0; j < n_u; ++j) {
                u_e += phi[j][qp] * c.get_elem_solution(_u_var)(j);
                v_e += phi[j][qp] * c.get_elem_solution(_v_var)(j);
            }
            Number u_nb = 0, v_nb = 0;
            for (std::size_t j = 0; j < n_u_n; ++j) {
                u_nb += phi_n[j][qp] * u_n_sol[j];
                v_nb += phi_n[j][qp] * v_n_sol[j];
            }

            // Pressures
            Number p_e = 0, p_nb = 0;
            for (std::size_t j = 0; j < n_p; ++j)
                p_e += psi[j][qp] * c.get_elem_solution(_p_var)(j);
            for (std::size_t j = 0; j < n_p_n; ++j)
                p_nb += psi_n[j][qp] * p_n_sol[j];

            // Average normal velocity: {u.n}
            const Number avg_un = 0.5 * ((u_e + u_nb) * n(0) + (v_e + v_nb) * n(1));
            const Number jump_p = p_e - p_nb;

            // Elem pressure test: -{q_e} * {u.n} + gamma_p * h * [[p]] * q_e
            for (std::size_t i = 0; i < n_p; ++i) {
                const Real psi_i = psi[i][qp];
                Fp(i) -= jxw * psi_i * avg_un;
                Fp(i) += jxw * Params::GAMMA_P * h_F * jump_p * psi_i;

                if (request_jacobian) {
                    for (std::size_t j = 0; j < n_u; ++j) {
                        (*Kpu_ee)(i,j) -= jxw * psi_i * 0.5 * phi[j][qp] * n(0);
                        (*Kpv_ee)(i,j) -= jxw * psi_i * 0.5 * phi[j][qp] * n(1);
                    }
                    for (std::size_t j = 0; j < n_u_n; ++j) {
                        (*Kpu_en)(i,j) -= jxw * psi_i * 0.5 * phi_n[j][qp] * n(0);
                        (*Kpv_en)(i,j) -= jxw * psi_i * 0.5 * phi_n[j][qp] * n(1);
                    }
                    for (std::size_t j = 0; j < n_p; ++j)
                        (*Kpp_ee)(i,j) += jxw * Params::GAMMA_P * h_F * psi[j][qp] * psi_i;
                    for (std::size_t j = 0; j < n_p_n; ++j)
                        (*Kpp_en)(i,j) -= jxw * Params::GAMMA_P * h_F * psi_n[j][qp] * psi_i;
                }
            }

            // Neighbor pressure test: +{q_n} * {u.n} - gamma_p * h * [[p]] * q_n
            for (std::size_t i = 0; i < n_p_n; ++i) {
                const Real psi_ni = psi_n[i][qp];
                Fp_n(i) += jxw * psi_ni * avg_un;
                Fp_n(i) -= jxw * Params::GAMMA_P * h_F * jump_p * psi_ni;

                if (request_jacobian) {
                    for (std::size_t j = 0; j < n_u; ++j) {
                        (*Kpu_ne)(i,j) += jxw * psi_ni * 0.5 * phi[j][qp] * n(0);
                        (*Kpv_ne)(i,j) += jxw * psi_ni * 0.5 * phi[j][qp] * n(1);
                    }
                    for (std::size_t j = 0; j < n_u_n; ++j) {
                        (*Kpu_nn)(i,j) += jxw * psi_ni * 0.5 * phi_n[j][qp] * n(0);
                        (*Kpv_nn)(i,j) += jxw * psi_ni * 0.5 * phi_n[j][qp] * n(1);
                    }
                    for (std::size_t j = 0; j < n_p; ++j)
                        (*Kpp_ne)(i,j) -= jxw * Params::GAMMA_P * h_F * psi[j][qp] * psi_ni;
                    for (std::size_t j = 0; j < n_p_n; ++j)
                        (*Kpp_nn)(i,j) += jxw * Params::GAMMA_P * h_F * psi_n[j][qp] * psi_ni;
                }
            }
        }
    } else {
        // ─── Boundary face: continuity term ──────────────────────────────────

        const BoundaryInfo& binfo = this->get_mesh().get_boundary_info();
        const auto bid_s = static_cast<unsigned short>(side);

        int bid_int = -1;
        for (int b : {Params::BID_INLET, Params::BID_OUTLET,
                      Params::BID_WALLS, Params::BID_CYLINDER})
            if (binfo.has_boundary_id(&elem, bid_s, static_cast<boundary_id_type>(b)))
                { bid_int = b; break; }
        if (bid_int < 0) return request_jacobian;

        if (bid_int == Params::BID_OUTLET) return request_jacobian;

        auto& Fp = c.get_elem_residual(_p_var);

        for (unsigned int qp = 0; qp < n_qp; ++qp) {
            const Real jxw = JxW[qp];
            const Point& n = normal[qp];
            const Point& xyz = u_fe->get_xyz()[qp];

            double u_bc = 0.0, v_bc = 0.0;
            if (bid_int == Params::BID_INLET)
                u_bc = inlet_u(xyz(1), this->time);

            const Number un_bc = u_bc * n(0) + v_bc * n(1);

            for (std::size_t i = 0; i < n_p; ++i)
                Fp(i) -= jxw * psi[i][qp] * un_bc;
        }
    }
    return request_jacobian;
}
