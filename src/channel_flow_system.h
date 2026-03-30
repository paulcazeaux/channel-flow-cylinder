#pragma once
/**
 * @file channel_flow_system.h
 * @brief FEMSystem subclass for steady 2D incompressible Navier-Stokes
 *        on the Schafer-Turek channel-cylinder domain (DFG 2D-1).
 *
 * Variables:
 *   u, v  — velocity components, P2 Lagrange (SECOND order)
 *   p     — pressure, P1 Lagrange (FIRST order)
 *
 * Weak form:
 *   Momentum:   ν∫∇u·∇w dx + ∫(u·∇)u·w dx − ∫p∇·w dx = 0
 *   Continuity: ∫q∇·u dx = 0
 *
 * Boundary conditions:
 *   Inlet   (BID 1): parabolic u = 4·U_MAX·y·(H−y)/H², v = 0
 *   Outlet  (BID 2): stress-free do-nothing (natural Neumann)
 *   Walls   (BID 3): no-slip u = v = 0
 *   Cylinder(BID 4): no-slip u = v = 0
 */

#include "libmesh/fem_system.h"
#include "libmesh/function_base.h"
#include "libmesh/mesh_base.h"

#include <memory>
#include <string>

/**
 * @class ChannelFlowSystem
 * @brief Steady incompressible Navier-Stokes via Taylor-Hood P2/P1 elements.
 */
class ChannelFlowSystem : public libMesh::FEMSystem
{
public:
    /**
     * @brief Construct and register with an EquationSystems instance.
     * @param es     Parent equation systems object.
     * @param name   System name used by libMesh.
     * @param number Unique system index assigned by libMesh.
     */
    ChannelFlowSystem(libMesh::EquationSystems& es,
                      const std::string& name,
                      unsigned int number);

    /**
     * @brief Add u, v (P2) and p (P1) variables; register Dirichlet BCs.
     *        Must be called before es.init().
     */
    void init_data() override;

    /**
     * @brief Request the FE data (JxW, phi, dphi) that element assembly needs.
     *        Called once per element before element_time_derivative /
     *        element_constraint.
     */
    void init_context(libMesh::DiffContext& ctx) override;

    /**
     * @brief Assemble the element momentum residual (and Jacobian if requested).
     *
     * Contributes ν∫∇u·∇w + (u·∇)u·w − p∇·w to the element residual.
     * Newton Jacobian is assembled when @p request_jacobian is true.
     *
     * @return Whether the Jacobian was assembled.
     */
    bool element_time_derivative(bool request_jacobian,
                                 libMesh::DiffContext& ctx) override;

    /**
     * @brief Assemble the element continuity residual (and Jacobian if requested).
     *
     * Contributes −∫q∇·u dx to the element residual.
     *
     * @return Whether the Jacobian was assembled.
     */
    bool element_constraint(bool request_jacobian,
                            libMesh::DiffContext& ctx) override;

    /// Variable index for x-velocity (set during init_data).
    unsigned int u_var() const { return _u_var; }
    /// Variable index for y-velocity (set during init_data).
    unsigned int v_var() const { return _v_var; }
    /// Variable index for pressure (set during init_data).
    unsigned int p_var() const { return _p_var; }

    /**
     * @brief Tag the outlet-corner node with BID_PRESSURE_PIN so that
     *        init_data() can apply a p=0 Dirichlet condition there.
     *
     * Must be called after mesh.all_second_order() and before
     * EquationSystems::init().  Removes the pressure null space of the
     * incompressible Stokes/NS system, which otherwise causes iterative
     * solvers to stall on the full saddle-point matrix.
     *
     * @param mesh The mesh to annotate (modified in place).
     */
    static void tag_pressure_pin(libMesh::MeshBase& mesh);

    /**
     * @brief Enable or disable Stokes (linear) mode.
     *
     * When true, advection terms (u·∇)u are suppressed in assembly, reducing
     * Navier-Stokes to the linear Stokes problem.  Newton converges in exactly
     * one iteration, providing a cheap initial guess for the full NS solve.
     */
    void set_stokes_mode(bool stokes) { _stokes_mode = stokes; }
    bool stokes_mode() const { return _stokes_mode; }

private:
    unsigned int _u_var; ///< P2 x-velocity variable index
    unsigned int _v_var; ///< P2 y-velocity variable index
    unsigned int _p_var; ///< P1 pressure variable index

    bool _stokes_mode = false; ///< Suppress advection for Stokes initialisation

    /// Owned parabolic inlet profile; kept alive for DirichletBoundary lifetime.
    std::unique_ptr<libMesh::FunctionBase<libMesh::Number>> _inlet_u_func;
};
