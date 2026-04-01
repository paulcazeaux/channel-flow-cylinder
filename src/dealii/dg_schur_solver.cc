/**
 * @file dg_schur_solver.cc
 * @brief Pressure Schur complement solver with Cahouet-Chabard preconditioner.
 *
 * System: [ A  B ] [ u ] = [ f ]    A = M_V + γ·dt·ν·K_DG (SPD)
 *         [ B^T 0] [ p ]   [ 0 ]    B = strong-form div (block-diagonal)
 *
 * Velocity solves use CG + BoomerAMG (strong threshold 0.1, Guermond).
 * Cahouet-Chabard: P_CC = ν M_p^{-1} + 1/(γdt) L_p^{-1} (both block-diagonal).
 */

#include "dg_navier_stokes.h"
#include "dg_params.h"

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_solver.h>

#include <iostream>

using namespace dealii;

// ═══════════════════════════════════════════════════════════════════════════════
//  AMG wrapper: converts deal.II SparseMatrix to PETSc and applies BoomerAMG
// ═══════════════════════════════════════════════════════════════════════════════

class AMGPreconditioner
{
public:
    AMGPreconditioner(const SparseMatrix<double>& matrix,
                      const SparsityPattern& sp,
                      double strong_threshold)
    {
        // Copy deal.II sparse matrix to PETSc using per-row nnz
        const unsigned int n = matrix.m();
        std::vector<unsigned int> row_lengths(n);
        for (unsigned int row = 0; row < n; ++row)
            row_lengths[row] = sp.row_length(row);

        petsc_matrix.reinit(n, n, row_lengths);

        for (unsigned int row = 0; row < n; ++row)
            for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                petsc_matrix.set(it->row(), it->column(), it->value());
        petsc_matrix.compress(VectorOperation::insert);

        // Configure BoomerAMG
        PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
        data.symmetric_operator = true;
        data.strong_threshold = strong_threshold;
        amg.initialize(petsc_matrix, data);
    }

    /// @brief Apply one V-cycle: dst = AMG^{-1} src
    void vmult(Vector<double>& dst, const Vector<double>& src) const
    {
        // Copy to PETSc vectors
        PETScWrappers::MPI::Vector petsc_src(MPI_COMM_SELF, src.size(), src.size());
        PETScWrappers::MPI::Vector petsc_dst(MPI_COMM_SELF, dst.size(), dst.size());
        for (unsigned int i = 0; i < src.size(); ++i)
            petsc_src(i) = src(i);
        petsc_src.compress(VectorOperation::insert);

        amg.vmult(petsc_dst, petsc_src);

        for (unsigned int i = 0; i < dst.size(); ++i)
            dst(i) = petsc_dst(i);
    }

private:
    PETScWrappers::SparseMatrix petsc_matrix;
    PETScWrappers::PreconditionBoomerAMG amg;
};

// ═══════════════════════════════════════════════════════════════════════════════
//  SchurOperator: -S = B^T A^{-1} B (positive semi-definite)
// ═══════════════════════════════════════════════════════════════════════════════

class SchurOperator
{
public:
    SchurOperator(const SparseMatrix<double>& A,
                  const AMGPreconditioner& A_amg,
                  const SparseMatrix<double>& B,
                  const SparseMatrix<double>& Bt,
                  unsigned int n_scalar_vel)
        : matrix_A(A), A_precond(A_amg)
        , matrix_B(B), matrix_Bt(Bt)
        , n_scalar(n_scalar_vel)
    {}

    void vmult(Vector<double>& dst, const Vector<double>& src) const
    {
        const unsigned int n_vel = 2 * n_scalar;
        Vector<double> Bp(n_vel), Ainv_Bp(n_vel);

        matrix_B.vmult(Bp, src);

        for (unsigned int d = 0; d < 2; ++d) {
            Vector<double> rhs_d(n_scalar), sol_d(n_scalar);
            for (unsigned int i = 0; i < n_scalar; ++i)
                rhs_d(i) = Bp(d * n_scalar + i);

            SolverControl control(100, 1e-3 * rhs_d.l2_norm());
            SolverCG<Vector<double>> cg(control);
            try {
                cg.solve(matrix_A, sol_d, rhs_d, A_precond);
            } catch (SolverControl::NoConvergence&) {}

            for (unsigned int i = 0; i < n_scalar; ++i)
                Ainv_Bp(d * n_scalar + i) = sol_d(i);
        }

        matrix_Bt.vmult(dst, Ainv_Bp);
        dst(0) = src(0);  // pin
    }

    unsigned int m() const { return matrix_Bt.m(); }
    unsigned int n() const { return matrix_Bt.m(); }

private:
    const SparseMatrix<double>& matrix_A;
    const AMGPreconditioner&    A_precond;
    const SparseMatrix<double>& matrix_B;
    const SparseMatrix<double>& matrix_Bt;
    unsigned int                n_scalar;
};

// ═══════════════════════════════════════════════════════════════════════════════
//  CahouetChabard: P_CC = ν M_p^{-1} + 1/(γdt) L_p^{-1}
//  Both are block-diagonal (DG element-local).
// ═══════════════════════════════════════════════════════════════════════════════

class CahouetChabard
{
public:
    CahouetChabard(double nu, double gamma_dt,
                   const std::vector<FullMatrix<double>>& mp_inv,
                   const std::vector<FullMatrix<double>>& lp_inv,
                   const DoFHandler<2>& pressure_dof_handler)
        : nu_(nu), gamma_dt_inv_(1.0 / gamma_dt)
        , mp_inv_(mp_inv), lp_inv_(lp_inv)
        , dof_handler_(pressure_dof_handler)
    {}

    void vmult(Vector<double>& dst, const Vector<double>& src) const
    {
        const unsigned int dpc = dof_handler_.get_fe().n_dofs_per_cell();
        std::vector<types::global_dof_index> dof_indices(dpc);
        Vector<double> local_src(dpc), local_mp(dpc), local_lp(dpc);

        dst = 0;
        unsigned int cell_idx = 0;
        for (const auto& cell : dof_handler_.active_cell_iterators()) {
            cell->get_dof_indices(dof_indices);
            for (unsigned int i = 0; i < dpc; ++i)
                local_src(i) = src(dof_indices[i]);

            mp_inv_[cell_idx].vmult(local_mp, local_src);
            lp_inv_[cell_idx].vmult(local_lp, local_src);

            for (unsigned int i = 0; i < dpc; ++i)
                dst(dof_indices[i]) = nu_ * local_mp(i)
                                    + gamma_dt_inv_ * local_lp(i);
            ++cell_idx;
        }
    }

private:
    double nu_, gamma_dt_inv_;
    const std::vector<FullMatrix<double>>& mp_inv_;
    const std::vector<FullMatrix<double>>& lp_inv_;
    const DoFHandler<2>& dof_handler_;
};

// ═══════════════════════════════════════════════════════════════════════════════
//  solve_schur
// ═══════════════════════════════════════════════════════════════════════════════

void DGNavierStokes::solve_schur(const Vector<double>& rhs_u,
                                 const Vector<double>& rhs_p,
                                 Vector<double>& sol_u,
                                 Vector<double>& sol_p)
{
    TimerOutput::Scope t(timer, "solve_schur");

    const unsigned int n_scalar = dof_handler_velocity.n_dofs();
    const double gamma_dt = Params::IMEX_GAMMA * Params::DT;

    // ── AMG on A (strong threshold 0.1, Guermond) ────────────────────────
    AMGPreconditioner A_amg(matrix_A, sparsity_A, /*strong_threshold=*/0.1);

    // ── Step 1: velocity pre-solve ───────────────────────────────────────
    Vector<double> z(n_vel_dofs);
    for (unsigned int d = 0; d < 2; ++d) {
        Vector<double> rhs_d(n_scalar), sol_d(n_scalar);
        for (unsigned int i = 0; i < n_scalar; ++i)
            rhs_d(i) = rhs_u(d * n_scalar + i);

        SolverControl control(200, 1e-8 * rhs_d.l2_norm());
        SolverCG<Vector<double>> cg(control);
        cg.solve(matrix_A, sol_d, rhs_d, A_amg);
        std::cout << "  Pre-solve (d=" << d << "): "
                  << control.last_step() << " CG iters\n";

        for (unsigned int i = 0; i < n_scalar; ++i)
            z(d * n_scalar + i) = sol_d(i);
    }

    // ── Step 2: Schur RHS ────────────────────────────────────────────────
    Vector<double> r(n_pres_dofs);
    matrix_Bt.vmult(r, z);
    r -= rhs_p;
    r(0) = 0.0;

    std::cout << "  Schur RHS: " << r.l2_norm() << "\n";

    // ── Step 3: pressure Schur solve ─────────────────────────────────────
    {
        SchurOperator schur_op(matrix_A, A_amg, matrix_B, matrix_Bt, n_scalar);

        // Build element-local L_p^{-1} blocks
        std::vector<FullMatrix<double>> lp_inv_blocks;
        {
            const unsigned int vel_dpc = fe_velocity_scalar.n_dofs_per_cell();
            const unsigned int pres_dpc = fe_pressure.n_dofs_per_cell();
            const QGauss<2> quad(Params::DG_ORDER + 1);
            FEValues<2> vel_fe(mapping, fe_velocity_scalar, quad,
                               update_gradients | update_JxW_values);
            FEValues<2> pres_fe(mapping, fe_pressure, quad,
                                update_values | update_JxW_values);

            lp_inv_blocks.resize(triangulation.n_active_cells());
            FullMatrix<double> local_B(vel_dpc, pres_dpc);
            FullMatrix<double> MvinvB(vel_dpc, pres_dpc);
            FullMatrix<double> local_Lp(pres_dpc, pres_dpc);

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
        }

        CahouetChabard precond(Params::NU, gamma_dt,
                               mp_inv_blocks, lp_inv_blocks,
                               dof_handler_pressure);

        SolverControl control(Params::SCHUR_MAX_IT,
                              1e-6 * r.l2_norm());
        SolverFGMRES<Vector<double>> fgmres(control);
        sol_p = 0;
        fgmres.solve(schur_op, sol_p, r, precond);
        sol_p(0) = 0.0;
        std::cout << "  Pressure: " << control.last_step() << " FGMRES iters\n";
    }

    // ── Step 4: velocity back-solve ──────────────────────────────────────
    {
        Vector<double> Bp(n_vel_dofs);
        matrix_B.vmult(Bp, sol_p);

        sol_u.reinit(n_vel_dofs);
        for (unsigned int d = 0; d < 2; ++d) {
            Vector<double> rhs_d(n_scalar), sol_d(n_scalar);
            for (unsigned int i = 0; i < n_scalar; ++i)
                rhs_d(i) = rhs_u(d * n_scalar + i) - Bp(d * n_scalar + i);

            SolverControl control(200, 1e-8 * rhs_d.l2_norm());
            SolverCG<Vector<double>> cg(control);
            cg.solve(matrix_A, sol_d, rhs_d, A_amg);
            std::cout << "  Back-solve (d=" << d << "): "
                      << control.last_step() << " CG iters\n";

            for (unsigned int i = 0; i < n_scalar; ++i)
                sol_u(d * n_scalar + i) = sol_d(i);
        }
    }
}
