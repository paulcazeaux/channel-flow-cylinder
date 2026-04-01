#pragma once
/**
 * @file channel_flow_system.h
 * @brief FEMSystem subclass for 2D incompressible Navier-Stokes using DG.
 *
 * Variables (all discontinuous, equal-order P_k):
 *   u, v  — velocity components, L2_LAGRANGE order DG_ORDER
 *   p     — pressure, L2_LAGRANGE order DG_ORDER
 *
 * Weak form (DG with SIPG viscous flux + upwind advection):
 *   Volume:   standard Galerkin integrals (same as CG)
 *   Faces:    SIPG penalty, pressure coupling, Lax-Friedrichs advection
 *   BCs:      weakly enforced via boundary face fluxes
 *   Pressure: equal-order stabilised with pressure-jump penalty
 */

#include "libmesh/fem_system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/newton_solver.h"

#include <memory>
#include <string>

/**
 * @class ChannelFlowSystem
 * @brief Incompressible Navier-Stokes via DG equal-order elements.
 */
class ChannelFlowSystem : public libMesh::FEMSystem
{
public:
    ChannelFlowSystem(libMesh::EquationSystems& es,
                      const std::string& name,
                      unsigned int number);

    /** @brief Add u, v, p as L2_LAGRANGE variables; enable internal sides. */
    void init_data() override;

    /** @brief Request FE data for element and side assembly. */
    void init_context(libMesh::DiffContext& ctx) override;

    /** @brief Build a DGFEMContext instead of the default FEMContext. */
    std::unique_ptr<libMesh::DiffContext> build_context() override;

    // ── Volume integrals (element interior) ──────────────────────────────────

    /** @brief Momentum residual: viscous + advection + pressure gradient. */
    bool element_time_derivative(bool request_jacobian,
                                 libMesh::DiffContext& ctx) override;

    /** @brief Continuity residual: -div(u). */
    bool element_constraint(bool request_jacobian,
                            libMesh::DiffContext& ctx) override;

    /** @brief Mass residual for time-dependent solves. */
    bool mass_residual(bool request_jacobian,
                       libMesh::DiffContext& ctx) override;

    // ── Face integrals (DG fluxes) ───────────────────────────────────────────

    /** @brief Side momentum fluxes: SIPG, upwind advection, pressure, weak BCs. */
    bool side_time_derivative(bool request_jacobian,
                              libMesh::DiffContext& ctx) override;

    /** @brief Side continuity fluxes: velocity divergence coupling + pressure jump. */
    bool side_constraint(bool request_jacobian,
                         libMesh::DiffContext& ctx) override;

    // ── Accessors ────────────────────────────────────────────────────────────

    unsigned int u_var() const { return _u_var; }
    unsigned int v_var() const { return _v_var; }
    unsigned int p_var() const { return _p_var; }

    /**
     * @brief Configure PETSc fieldsplit preconditioner for DG DOFs.
     *
     * DG DOFs are element-local (no sharing); collects velocity and pressure
     * index sets by iterating over elements.
     */
    static void configure_fieldsplit(ChannelFlowSystem& sys,
                                     libMesh::MeshBase& mesh,
                                     libMesh::NewtonSolver& newton);

    /**
     * @brief Set up PETSc MatNullSpace for pressure constant mode.
     *
     * Must be called after the first assembly so the matrix exists.
     * Projects out the constant pressure mode to remove the null space.
     */
    static void setup_pressure_null_space(ChannelFlowSystem& sys,
                                          libMesh::MeshBase& mesh,
                                          libMesh::NewtonSolver& newton);

    void set_stokes_mode(bool stokes) { _stokes_mode = stokes; }
    bool stokes_mode() const { return _stokes_mode; }

private:
    unsigned int _u_var;
    unsigned int _v_var;
    unsigned int _p_var;

    bool _stokes_mode = false;
};
