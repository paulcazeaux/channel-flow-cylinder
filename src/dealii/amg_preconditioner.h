#pragma once
/**
 * @file amg_preconditioner.h
 * @brief AMG wrapper: converts deal.II SparseMatrix to PETSc and applies BoomerAMG.
 */

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

class AMGPreconditioner
{
public:
    AMGPreconditioner(const dealii::SparseMatrix<double>& matrix,
                      const dealii::SparsityPattern& sp,
                      double strong_threshold)
    {
        const unsigned int n = matrix.m();
        std::vector<unsigned int> row_lengths(n);
        for (unsigned int row = 0; row < n; ++row)
            row_lengths[row] = sp.row_length(row);

        petsc_matrix.reinit(n, n, row_lengths);

        for (unsigned int row = 0; row < n; ++row)
            for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                petsc_matrix.set(it->row(), it->column(), it->value());
        petsc_matrix.compress(dealii::VectorOperation::insert);

        dealii::PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
        data.symmetric_operator = true;
        data.strong_threshold = strong_threshold;
        amg.initialize(petsc_matrix, data);

        n_dofs = n;
    }

    void vmult(dealii::Vector<double>& dst,
               const dealii::Vector<double>& src) const
    {
        dealii::PETScWrappers::MPI::Vector ps(MPI_COMM_SELF, n_dofs, n_dofs);
        dealii::PETScWrappers::MPI::Vector pd(MPI_COMM_SELF, n_dofs, n_dofs);
        for (unsigned int i = 0; i < n_dofs; ++i)
            ps(i) = src(i);
        ps.compress(dealii::VectorOperation::insert);

        amg.vmult(pd, ps);

        for (unsigned int i = 0; i < n_dofs; ++i)
            dst(i) = pd(i);
    }

private:
    dealii::PETScWrappers::SparseMatrix petsc_matrix;
    dealii::PETScWrappers::PreconditionBoomerAMG amg;
    unsigned int n_dofs;
};
