/**
 * @file drag_lift.cpp
 * @brief Drag and lift force computation by boundary integration over BID_CYLINDER.
 *
 * Integrates the viscous stress traction σ·n over all sides of active elements
 * that lie on the cylinder boundary.  Uses side FE reinit so that shape functions
 * and normals are evaluated at the correct quadrature points on each edge.
 */

#include "drag_lift.h"
#include "channel_flow_system.h"
#include "params.h"

#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel.h"

#include <vector>

std::pair<double, double> compute_drag_lift(const ChannelFlowSystem& sys,
                                            const libMesh::MeshBase& mesh)
{
    const libMesh::DofMap& dof_map = sys.get_dof_map();

    // ── Build side FE objects ─────────────────────────────────────────────────
    // QGauss of dimension 1 integrates over edges in 2D.
    // Attach to 2D FEBase; reinit(elem, side) projects to the correct edge.
    const libMesh::FEType vel_type = sys.variable_type(sys.u_var());
    const libMesh::FEType p_type   = sys.variable_type(sys.p_var());

    auto vel_fe = libMesh::FEBase::build(2, vel_type);
    auto p_fe   = libMesh::FEBase::build(2, p_type);

    libMesh::QGauss qrule(1, vel_type.default_quadrature_order());
    vel_fe->attach_quadrature_rule(&qrule);
    p_fe->attach_quadrature_rule(&qrule);

    const auto& JxW    = vel_fe->get_JxW();
    const auto& phi    = vel_fe->get_phi();
    const auto& dphi   = vel_fe->get_dphi();
    const auto& normal = vel_fe->get_normals();
    const auto& p_phi  = p_fe->get_phi();

    // ── Solution vector ───────────────────────────────────────────────────────
    // current_local_solution is populated (with ghost entries) after sys.solve().
    const libMesh::NumericVector<libMesh::Number>& sol =
        *sys.current_local_solution;

    const double nu = Params::NU;
    const libMesh::boundary_id_type cyl_bid =
        static_cast<libMesh::boundary_id_type>(Params::BID_CYLINDER);

    // ── Boundary integration ──────────────────────────────────────────────────
    double F_D = 0.0, F_L = 0.0;

    std::vector<libMesh::dof_id_type> u_dofs, v_dofs, p_dofs;

    for (const auto* elem : mesh.active_local_element_ptr_range()) {
        for (unsigned int s = 0; s < elem->n_sides(); ++s) {
            if (!mesh.get_boundary_info().has_boundary_id(elem, s, cyl_bid))
                continue;

            vel_fe->reinit(elem, s);
            p_fe->reinit(elem, s);

            dof_map.dof_indices(elem, u_dofs, sys.u_var());
            dof_map.dof_indices(elem, v_dofs, sys.v_var());
            dof_map.dof_indices(elem, p_dofs, sys.p_var());

            const unsigned int n_qp = qrule.n_points();

            for (unsigned int qp = 0; qp < n_qp; ++qp) {
                // Interpolate p, ∂u/∂x, ∂u/∂y, ∂v/∂x, ∂v/∂y at quadrature point
                double p_qp   = 0.0;
                double du_dx  = 0.0, du_dy = 0.0;
                double dv_dx  = 0.0, dv_dy = 0.0;

                for (std::size_t j = 0; j < p_dofs.size(); ++j)
                    p_qp += p_phi[j][qp] * sol(p_dofs[j]);

                for (std::size_t j = 0; j < u_dofs.size(); ++j) {
                    const double u_j = sol(u_dofs[j]);
                    const double v_j = sol(v_dofs[j]);
                    du_dx += dphi[j][qp](0) * u_j;
                    du_dy += dphi[j][qp](1) * u_j;
                    dv_dx += dphi[j][qp](0) * v_j;
                    dv_dy += dphi[j][qp](1) * v_j;
                }

                // n_libMesh is the outward normal of the fluid element (fluid→cylinder).
                // The Schafer-Turek convention uses n_body = -n_libMesh (cylinder→fluid).
                // F = ∫ σ · n_body dS = -∫ σ · n_libMesh dS
                const libMesh::Point& n = normal[qp];

                // Stress traction σ·n_libMesh  (σ = -pI + ν(∇u + ∇uᵀ))
                const double tx = -p_qp * n(0)
                                  + nu * (2.0 * du_dx * n(0)
                                          + (du_dy + dv_dx) * n(1));
                const double ty = -p_qp * n(1)
                                  + nu * ((du_dy + dv_dx) * n(0)
                                          + 2.0 * dv_dy * n(1));

                F_D += JxW[qp] * tx;
                F_L += JxW[qp] * ty;
            }
        }
    }

    // Reduce across MPI ranks
    sys.comm().sum(F_D);
    sys.comm().sum(F_L);

    // Negate: libMesh normals point fluid→cylinder; Schafer-Turek convention
    // uses the outward body normal (cylinder→fluid), so F = -∫ σ·n_libMesh dS.
    return {-F_D, -F_L};
}
