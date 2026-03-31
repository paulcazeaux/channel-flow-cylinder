/**
 * @file channel_flow_assembly.cpp
 * @brief Element-level weak-form assembly for ChannelFlowSystem.
 *
 * Implements the Taylor-Hood P2/P1 Galerkin form for incompressible
 * Navier-Stokes, including the mass matrix for time-dependent solves.
 * Variable setup and Dirichlet BCs are in channel_flow_system.cpp.
 */

#include "channel_flow_system.h"
#include "params.h"

#include "libmesh/fem_context.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature.h"
#include "libmesh/dof_map.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// ── Element momentum assembly ─────────────────────────────────────────────────

bool ChannelFlowSystem::element_time_derivative(bool request_jacobian,
                                                libMesh::DiffContext& ctx)
{
    libMesh::FEMContext& c = libMesh::cast_ref<libMesh::FEMContext&>(ctx);

    libMesh::FEBase* u_fe = nullptr;
    libMesh::FEBase* p_fe = nullptr;
    c.get_element_fe(_u_var, u_fe);
    c.get_element_fe(_p_var, p_fe);

    const std::vector<libMesh::Real>&                          JxW  = u_fe->get_JxW();
    const std::vector<std::vector<libMesh::Real>>&             phi  = u_fe->get_phi();
    const std::vector<std::vector<libMesh::RealGradient>>&     dphi = u_fe->get_dphi();
    const std::vector<std::vector<libMesh::Real>>&             psi  = p_fe->get_phi();

    const std::size_t n_u_dofs = c.get_dof_indices(_u_var).size();
    const std::size_t n_p_dofs = c.get_dof_indices(_p_var).size();
    const unsigned int n_qp    = c.get_element_qrule().n_points();

    libMesh::DenseSubVector<libMesh::Number>& Fu = c.get_elem_residual(_u_var);
    libMesh::DenseSubVector<libMesh::Number>& Fv = c.get_elem_residual(_v_var);

    libMesh::DenseSubMatrix<libMesh::Number>* Kuu = nullptr;
    libMesh::DenseSubMatrix<libMesh::Number>* Kuv = nullptr;
    libMesh::DenseSubMatrix<libMesh::Number>* Kup = nullptr;
    libMesh::DenseSubMatrix<libMesh::Number>* Kvu = nullptr;
    libMesh::DenseSubMatrix<libMesh::Number>* Kvv = nullptr;
    libMesh::DenseSubMatrix<libMesh::Number>* Kvp = nullptr;
    if (request_jacobian) {
        Kuu = &c.get_elem_jacobian(_u_var, _u_var);
        Kuv = &c.get_elem_jacobian(_u_var, _v_var);
        Kup = &c.get_elem_jacobian(_u_var, _p_var);
        Kvu = &c.get_elem_jacobian(_v_var, _u_var);
        Kvv = &c.get_elem_jacobian(_v_var, _v_var);
        Kvp = &c.get_elem_jacobian(_v_var, _p_var);
    }

    const double nu = Params::NU;

    for (unsigned int qp = 0; qp < n_qp; ++qp) {
        const libMesh::Number       u      = c.interior_value(_u_var, qp);
        const libMesh::Number       v      = c.interior_value(_v_var, qp);
        const libMesh::Number       p      = c.interior_value(_p_var, qp);
        const libMesh::Gradient     grad_u = c.interior_gradient(_u_var, qp);
        const libMesh::Gradient     grad_v = c.interior_gradient(_v_var, qp);
        const libMesh::Real         jxw    = JxW[qp];

        for (std::size_t i = 0; i < n_u_dofs; ++i) {
            const libMesh::Real          phi_i  = phi[i][qp];
            const libMesh::RealGradient& dphi_i = dphi[i][qp];

            // Momentum (x): ν ∇u·∇w − p ∂w/∂x  [+ advection if !stokes_mode]
            Fu(i) += jxw * (nu * (grad_u * dphi_i) - p * dphi_i(0));
            Fv(i) += jxw * (nu * (grad_v * dphi_i) - p * dphi_i(1));
            if (!_stokes_mode) {
                Fu(i) += jxw * (u * grad_u(0) + v * grad_u(1)) * phi_i;
                Fv(i) += jxw * (u * grad_v(0) + v * grad_v(1)) * phi_i;
            }

            if (request_jacobian) {
                for (std::size_t j = 0; j < n_u_dofs; ++j) {
                    const libMesh::Real          phi_j  = phi[j][qp];
                    const libMesh::RealGradient& dphi_j = dphi[j][qp];

                    // Viscous contributions (always present)
                    (*Kuu)(i,j) += jxw * nu * (dphi_j * dphi_i);
                    (*Kvv)(i,j) += jxw * nu * (dphi_j * dphi_i);

                    // Advection Jacobian (suppressed in Stokes mode)
                    if (!_stokes_mode) {
                        (*Kuu)(i,j) += jxw * (phi_j * grad_u(0) + u * dphi_j(0) + v * dphi_j(1)) * phi_i;
                        (*Kuv)(i,j) += jxw * phi_j * grad_u(1) * phi_i;
                        (*Kvu)(i,j) += jxw * phi_j * grad_v(0) * phi_i;
                        (*Kvv)(i,j) += jxw * (u * dphi_j(0) + phi_j * grad_v(1) + v * dphi_j(1)) * phi_i;
                    }
                }
                for (std::size_t j = 0; j < n_p_dofs; ++j) {
                    const libMesh::Real psi_j = psi[j][qp];
                    (*Kup)(i,j) += jxw * (-psi_j * dphi_i(0));
                    (*Kvp)(i,j) += jxw * (-psi_j * dphi_i(1));
                }
            }
        }
    }
    return request_jacobian;
}

// ── Element continuity assembly ───────────────────────────────────────────────

bool ChannelFlowSystem::element_constraint(bool request_jacobian,
                                           libMesh::DiffContext& ctx)
{
    libMesh::FEMContext& c = libMesh::cast_ref<libMesh::FEMContext&>(ctx);

    libMesh::FEBase* u_fe = nullptr;
    libMesh::FEBase* p_fe = nullptr;
    c.get_element_fe(_u_var, u_fe);
    c.get_element_fe(_p_var, p_fe);

    const std::vector<libMesh::Real>&                      JxW  = u_fe->get_JxW();
    const std::vector<std::vector<libMesh::RealGradient>>& dphi = u_fe->get_dphi();
    const std::vector<std::vector<libMesh::Real>>&         psi  = p_fe->get_phi();

    const std::size_t n_u_dofs = c.get_dof_indices(_u_var).size();
    const std::size_t n_p_dofs = c.get_dof_indices(_p_var).size();
    const unsigned int n_qp    = c.get_element_qrule().n_points();

    libMesh::DenseSubVector<libMesh::Number>& Fp = c.get_elem_residual(_p_var);

    libMesh::DenseSubMatrix<libMesh::Number>* Kpu = nullptr;
    libMesh::DenseSubMatrix<libMesh::Number>* Kpv = nullptr;
    if (request_jacobian) {
        Kpu = &c.get_elem_jacobian(_p_var, _u_var);
        Kpv = &c.get_elem_jacobian(_p_var, _v_var);
    }

    for (unsigned int qp = 0; qp < n_qp; ++qp) {
        const libMesh::Gradient grad_u = c.interior_gradient(_u_var, qp);
        const libMesh::Gradient grad_v = c.interior_gradient(_v_var, qp);
        const libMesh::Real     jxw   = JxW[qp];

        for (std::size_t i = 0; i < n_p_dofs; ++i) {
            const libMesh::Real psi_i = psi[i][qp];

            // Continuity: −∫ q ∇·u dx
            Fp(i) -= jxw * (grad_u(0) + grad_v(1)) * psi_i;

            if (request_jacobian) {
                for (std::size_t j = 0; j < n_u_dofs; ++j) {
                    (*Kpu)(i,j) -= jxw * dphi[j][qp](0) * psi_i;
                    (*Kpv)(i,j) -= jxw * dphi[j][qp](1) * psi_i;
                }
            }
        }
    }
    return request_jacobian;
}

// ── Mass matrix assembly (time-dependent) ────────────────────────────────────

bool ChannelFlowSystem::mass_residual(bool request_jacobian,
                                      libMesh::DiffContext& ctx)
{
    libMesh::FEMContext& c = libMesh::cast_ref<libMesh::FEMContext&>(ctx);

    libMesh::FEBase* u_fe = nullptr;
    c.get_element_fe(_u_var, u_fe);

    const std::vector<libMesh::Real>&              JxW = u_fe->get_JxW();
    const std::vector<std::vector<libMesh::Real>>& phi = u_fe->get_phi();

    const std::size_t n_u_dofs = c.get_dof_indices(_u_var).size();
    const unsigned int n_qp    = c.get_element_qrule().n_points();

    libMesh::DenseSubVector<libMesh::Number>& Fu = c.get_elem_residual(_u_var);
    libMesh::DenseSubVector<libMesh::Number>& Fv = c.get_elem_residual(_v_var);

    libMesh::DenseSubMatrix<libMesh::Number>* Muu = nullptr;
    libMesh::DenseSubMatrix<libMesh::Number>* Mvv = nullptr;
    if (request_jacobian) {
        Muu = &c.get_elem_jacobian(_u_var, _u_var);
        Mvv = &c.get_elem_jacobian(_v_var, _v_var);
    }

    for (unsigned int qp = 0; qp < n_qp; ++qp) {
        const libMesh::Number u_dot = c.interior_value(_u_var, qp);
        const libMesh::Number v_dot = c.interior_value(_v_var, qp);
        const libMesh::Real   jxw   = JxW[qp];

        for (std::size_t i = 0; i < n_u_dofs; ++i) {
            const libMesh::Real phi_i = phi[i][qp];

            // ∫ ρ u̇ · w dx  (ρ = 1 for incompressible)
            Fu(i) += jxw * u_dot * phi_i;
            Fv(i) += jxw * v_dot * phi_i;

            if (request_jacobian) {
                for (std::size_t j = 0; j < n_u_dofs; ++j) {
                    const libMesh::Real phi_j = phi[j][qp];
                    (*Muu)(i,j) += jxw * phi_j * phi_i;
                    (*Mvv)(i,j) += jxw * phi_j * phi_i;
                }
            }
        }
    }
    return request_jacobian;
}
