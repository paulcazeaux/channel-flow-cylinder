#pragma once
/**
 * @file petsc_utils.h
 * @brief Utility helpers for PETSc option strings.
 */

#include <iomanip>
#include <sstream>
#include <string>

/**
 * @brief Format a double for PetscOptionsSetValue using scientific notation.
 *
 * std::to_string uses %f (6 decimal places fixed) and silently truncates values
 * like 1e-10 to "0.000000".  PETSc then reads rtol=0, a criterion that is never
 * satisfied, causing the solver to always exhaust max_it and return DIVERGED_ITS.
 */
inline std::string petsc_str(double v)
{
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(6) << v;
    return ss.str();
}
