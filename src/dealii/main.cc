/**
 * @file main.cc
 * @brief Entry point for the DG Navier-Stokes solver (deal.II).
 */

#include "dg_navier_stokes.h"

#include <deal.II/base/utilities.h>

#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    dealii::Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

    DGNavierStokes solver;
    solver.run();

    return 0;
}
