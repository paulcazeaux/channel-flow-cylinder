/**
 * @file channel_flow_system.cpp
 * @brief ChannelFlowSystem: variable setup, Dirichlet BCs, and mesh utilities.
 *
 * Weak-form assembly is in channel_flow_assembly.cpp.
 */

#include "channel_flow_system.h"
#include "params.h"

#include "libmesh/dof_map.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/point.h"
#include "libmesh/dense_vector.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/boundary_info.h"
#include "libmesh/node.h"

#include <limits>
#include <set>

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

// ── Static helpers ────────────────────────────────────────────────────────────

void ChannelFlowSystem::tag_pressure_pin(libMesh::MeshBase& mesh)
{
    // Find the local node nearest to the outlet-bottom corner (CHANNEL_LENGTH, 0).
    // That node is tagged BID_PRESSURE_PIN so init_data() can add p=0 Dirichlet
    // there, removing the pressure null space for iterative solvers.
    libMesh::Node* pin_node = nullptr;
    double min_dist2 = std::numeric_limits<double>::max();
    for (auto& node : mesh.local_node_ptr_range()) {
        const double dx = (*node)(0) - Params::CHANNEL_LENGTH;
        const double dy = (*node)(1);
        if (dx * dx + dy * dy < min_dist2) {
            min_dist2 = dx * dx + dy * dy;
            pin_node  = node;
        }
    }
    libmesh_assert_msg(pin_node, "tag_pressure_pin: no local nodes found");
    mesh.get_boundary_info().add_node(
        pin_node,
        static_cast<libMesh::boundary_id_type>(Params::BID_PRESSURE_PIN));
}

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

    // Pressure pin: p=0 at BID_PRESSURE_PIN (outlet corner).
    // Fixes the pressure null space so iterative solvers can converge on the
    // saddle-point system.  Physically valid: incompressible pressure is only
    // determined up to a constant; pinning one node sets the reference level.
    {
        const std::set<libMesh::boundary_id_type> ids = {
            static_cast<libMesh::boundary_id_type>(Params::BID_PRESSURE_PIN)
        };
        libMesh::ZeroFunction<libMesh::Number> zero;
        this->get_dof_map().add_dirichlet_boundary(
            libMesh::DirichletBoundary(ids, {_p_var}, &zero));
    }

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
