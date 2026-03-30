/**
 * @file channel_flow_system.cpp
 * @brief Implementation of ChannelFlowSystem: variable setup, Dirichlet BCs,
 *        and Taylor-Hood weak-form assembly for steady Navier-Stokes.
 */

#include "channel_flow_system.h"
#include "params.h"

#include "libmesh/fem_context.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature.h"
#include "libmesh/dof_map.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/dense_vector.h"
#include "libmesh/point.h"
#include "libmesh/libmesh_common.h"

#include <set>
#include <vector>
#include <sstream>

// ── Inlet velocity function ───────────────────────────────────────────────────

namespace {

/**
 * @brief Parabolic inlet u-velocity: u(y) = 4·U_MAX·y·(H−y)/H²
 *
 * Implements FunctionBase so it can be passed to DirichletBoundary.
 * DirichletBoundary clones this object internally, so no lifetime concern.
 */
class InletVelocityU : public libMesh::FunctionBase<libMesh::Number>
{
public:
    std::unique_ptr<libMesh::FunctionBase<libMesh::Number>>
    clone() const override
    {
        return std::make_unique<InletVelocityU>();
    }

    libMesh::Number operator()(const libMesh::Point& p,
                               const libMesh::Real /*t*/ = 0.0) override
    {
        const double y = p(1);
        return 4.0 * Params::U_MAX * y
               * (Params::CHANNEL_HEIGHT - y)
               / (Params::CHANNEL_HEIGHT * Params::CHANNEL_HEIGHT);
    }

    void operator()(const libMesh::Point& p,
                    const libMesh::Real t,
                    libMesh::DenseVector<libMesh::Number>& output) override
    {
        output.resize(1);
        output(0) = (*this)(p, t);
    }
};

} // anonymous namespace

// ── ChannelFlowSystem ─────────────────────────────────────────────────────────

ChannelFlowSystem::ChannelFlowSystem(libMesh::EquationSystems& es,
                                     const std::string& name,
                                     unsigned int number)
    : libMesh::FEMSystem(es, name, number),
      _u_var(0), _v_var(0), _p_var(0)
{}

void ChannelFlowSystem::init_data()
{
    // ── Add Taylor-Hood variables ─────────────────────────────────────────────
    _u_var = this->add_variable("u", libMesh::SECOND, libMesh::LAGRANGE);
    _v_var = this->add_variable("v", libMesh::SECOND, libMesh::LAGRANGE);
    _p_var = this->add_variable("p", libMesh::FIRST,  libMesh::LAGRANGE);

    // ── Dirichlet BCs ─────────────────────────────────────────────────────────

    // No-slip: top/bottom walls and cylinder surface
    {
        const std::set<libMesh::boundary_id_type> ids = {
            static_cast<libMesh::boundary_id_type>(Params::BID_WALLS),
            static_cast<libMesh::boundary_id_type>(Params::BID_CYLINDER)
        };
        const std::vector<unsigned int> vel_vars = {_u_var, _v_var};
        libMesh::ZeroFunction<libMesh::Number> zero;
        this->get_dof_map().add_dirichlet_boundary(
            libMesh::DirichletBoundary(ids, vel_vars, &zero));
    }

    // Inlet u-velocity: parabolic profile
    {
        const std::set<libMesh::boundary_id_type> ids = {
            static_cast<libMesh::boundary_id_type>(Params::BID_INLET)
        };
        InletVelocityU inlet_u;
        this->get_dof_map().add_dirichlet_boundary(
            libMesh::DirichletBoundary(ids, {_u_var}, &inlet_u));
    }

    // Inlet v-velocity: zero
    {
        const std::set<libMesh::boundary_id_type> ids = {
            static_cast<libMesh::boundary_id_type>(Params::BID_INLET)
        };
        libMesh::ZeroFunction<libMesh::Number> zero;
        this->get_dof_map().add_dirichlet_boundary(
            libMesh::DirichletBoundary(ids, {_v_var}, &zero));
    }

    // Outlet (BID_OUTLET): natural Neumann — no action needed.

    // ── Delegate to parent to finalize variable/DOF setup ─────────────────────
    this->FEMSystem::init_data();
}

void ChannelFlowSystem::init_context(libMesh::DiffContext& ctx)
{
    libMesh::FEMContext& c = libMesh::cast_ref<libMesh::FEMContext&>(ctx);

    libMesh::FEBase* u_fe = nullptr;
    libMesh::FEBase* p_fe = nullptr;
    c.get_element_fe(_u_var, u_fe);
    c.get_element_fe(_p_var, p_fe);

    // Request the data arrays actually used in element assembly.
    u_fe->get_JxW();
    u_fe->get_phi();
    u_fe->get_dphi();
    p_fe->get_phi();
    // p_fe dphi not needed (continuity uses velocity shape gradients only).
}

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
