#pragma once
/**
 * @file dg_navier_stokes.h
 * @brief DG incompressible Navier-Stokes solver using deal.II.
 *
 * SIPG viscous + Lax-Friedrichs advection + equal-order Q_k/Q_k with
 * pressure-jump stabilization.  IMEX-ARK3 time integration.  Pressure
 * Schur complement solver with Cahouet-Chabard preconditioner.
 */

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <array>
#include <string>
#include <vector>

/**
 * @class DGNavierStokes
 * @brief 2D incompressible Navier-Stokes with DG and IMEX time stepping.
 *
 * Uses separate DoFHandlers for velocity (vector, 2 components) and
 * pressure (scalar) to enable clean block assembly and Schur complement
 * solves.  All matrices (A, B, M_V, M_p, L_p) are assembled separately.
 */
class DGNavierStokes
{
public:
    DGNavierStokes();

    /// @brief Run the full simulation.
    void run();

private:
    // ── Setup ────────────────────────────────────────────────────────────────

    /// @brief Generate quad mesh for channel with cylinder (Schafer-Turek).
    void setup_mesh(unsigned int n_refinements = 2);

    /// @brief Create DG finite element spaces and distribute DOFs.
    void setup_dofs();

    /// @brief Allocate and compute sparsity patterns for all block matrices.
    void setup_sparsity();

    /// @brief Pre-compute element-local mass matrix inverses (DG block-diagonal).
    void setup_mass_inverses();

    // ── Assembly ─────────────────────────────────────────────────────────────

    /// @brief Assemble the implicit operator A = M_V + gamma*dt*nu*K_DG.
    void assemble_implicit_operator();

    /// @brief Assemble A = nu*K_DG only (no mass term) for steady Stokes test.
    void assemble_viscous_only();

    /// @brief Assemble pressure gradient B and divergence B^T operators.
    void assemble_pressure_coupling();

    /// @brief Assemble pressure-jump stabilization C.
    void assemble_pressure_stabilization();

    /// @brief Assemble discrete pressure Laplacian L_p = B^T M_V^{-1} B.
    void assemble_pressure_laplacian();

    /// @brief Evaluate explicit advection: N(u) with Lax-Friedrichs flux.
    void evaluate_advection(const dealii::Vector<double>& velocity,
                            dealii::Vector<double>& advection_rhs,
                            double time = 0.0);

    // ── Solve ────────────────────────────────────────────────────────────────

    /// @brief Solve one Schur complement system for (u, p).
    void solve_schur(const dealii::Vector<double>& rhs_u,
                     const dealii::Vector<double>& rhs_p,
                     dealii::Vector<double>& sol_u,
                     dealii::Vector<double>& sol_p);

    // ── Output ───────────────────────────────────────────────────────────────

    /// @brief Assemble SIPG weak Dirichlet RHS for velocity.
    void assemble_bc_rhs(dealii::Vector<double>& rhs_u,
                         double coeff, double t) const;

    /// @brief Apply block-diagonal M_V to a velocity vector.
    void apply_mass(const dealii::Vector<double>& src,
                    dealii::Vector<double>& dst) const;

    /// @brief Write VTK output for the current time step.
    void output_results(double time, int step) const;

    /// @brief Compute drag and lift on the cylinder boundary.
    std::pair<double, double> compute_drag_lift() const;

    // ── Inlet BC ─────────────────────────────────────────────────────────────

    /// @brief Parabolic inlet u-velocity with smooth ramp.
    static double inlet_u(double y, double t);

    // ── Data members ─────────────────────────────────────────────────────────

    dealii::Triangulation<2> triangulation;

    dealii::MappingQ<2>   mapping;              ///< Higher-order mapping for curved cells
    dealii::FE_DGQ<2>     fe_velocity_scalar;   ///< Scalar DG basis for one velocity component (Q_k)
    dealii::FE_DGQ<2>     fe_pressure;           ///< Scalar DG basis for pressure (Q_{k-1})
    dealii::DoFHandler<2> dof_handler_velocity; ///< DOFs for (u, v) — 2 × scalar
    dealii::DoFHandler<2> dof_handler_pressure; ///< DOFs for p

    unsigned int n_vel_dofs;   ///< Total velocity DOFs (2 × n_scalar_vel_dofs)
    unsigned int n_pres_dofs;  ///< Total pressure DOFs

    // Block matrices (deal.II native sparse)
    dealii::SparsityPattern      sparsity_A, sparsity_B, sparsity_Bt, sparsity_C, sparsity_Lp;
    dealii::SparseMatrix<double> matrix_A;   ///< M_V + gamma*dt*nu*K_DG (SPD)
    dealii::SparseMatrix<double> matrix_B;   ///< Pressure gradient (n_vel × n_pres)
    dealii::SparseMatrix<double> matrix_Bt;  ///< Velocity divergence (n_pres × n_vel)
    dealii::SparseMatrix<double> matrix_C;   ///< Pressure-jump stabilization
    dealii::SparseMatrix<double> matrix_Lp;  ///< B^T M_V^{-1} B (pressure Laplacian)
    dealii::SparseMatrix<double> mass_velocity; ///< M_V (velocity mass)
    dealii::SparseMatrix<double> mass_pressure; ///< M_p (pressure mass)

    /// Element-local velocity mass inverse blocks: M_V_inv[elem] is dense (n_local × n_local).
    std::vector<dealii::FullMatrix<double>> mv_inv_blocks;

    /// Element-local pressure mass inverse blocks.
    std::vector<dealii::FullMatrix<double>> mp_inv_blocks;

    // Solution vectors
    dealii::Vector<double> solution_velocity;  ///< Current velocity (u, v interleaved)
    dealii::Vector<double> solution_pressure;  ///< Current pressure
    dealii::Vector<double> old_velocity;       ///< Previous time step velocity

    /// PVD time series entries (accumulated across output calls).
    mutable std::vector<std::pair<double, std::string>> pvd_entries;

    dealii::TimerOutput timer;
};
