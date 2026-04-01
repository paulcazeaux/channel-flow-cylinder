#pragma once
/**
 * @file dg_params.h
 * @brief Centralized parameters for the DG channel-flow solver (deal.II).
 *
 * Physical constants match the libMesh solver (Schafer-Turek DFG 2D-1/2D-2).
 * DG-specific parameters (penalty, stabilization, IMEX) are added here.
 */

namespace Params {

// ── Physical domain (Schafer-Turek 2D-1) ─────────────────────────────────────
constexpr double CHANNEL_LENGTH = 2.2;
constexpr double CHANNEL_HEIGHT = 0.41;
constexpr double CYL_X          = 0.2;
constexpr double CYL_Y          = 0.2;
constexpr double CYL_RADIUS     = 0.05;
constexpr double CYL_DIAMETER   = 2.0 * CYL_RADIUS;

// ── Fluid properties ──────────────────────────────────────────────────────────
constexpr double NU  = 1.0e-3;
constexpr double RHO = 1.0;

// ── Inlet velocity ───────────────────────────────────────────────────────────
constexpr double U_MAX  = 1.5;
constexpr double U_MEAN = (2.0 / 3.0) * U_MAX;
constexpr double RE     = U_MEAN * CYL_DIAMETER / NU;

// ── Boundary IDs (deal.II GridGenerator::channel_with_cylinder, colorize=true)
constexpr unsigned int BID_INLET    = 0;
constexpr unsigned int BID_OUTLET   = 1;
constexpr unsigned int BID_CYLINDER = 2;
constexpr unsigned int BID_WALLS    = 3;

// ── DG discretization ────────────────────────────────────────────────────────
constexpr unsigned int DG_ORDER = 2;  ///< Velocity order k; pressure uses k-1

/// SIPG viscous penalty parameter (dimensionless).
constexpr double SIGMA_V = 10.0 * DG_ORDER * DG_ORDER;

/// Pressure-jump stabilization (not needed for Q_k/Q_{k-1}, set to 0).
constexpr double GAMMA_P = 0.0;

// ── Time stepping (IMEX-ARK3) ────────────────────────────────────────────────
constexpr double DT      = 0.001;
constexpr double T_FINAL = 0.1;   // Short run for testing
constexpr double T_RAMP  = 0.05;  // Quick ramp for testing

/// SDIRK diagonal coefficient for IMEX-ARK3 (Kennedy & Carpenter 2003).
constexpr double IMEX_GAMMA = 0.4358665215;

// ── Linear solver ────────────────────────────────────────────────────────────
constexpr int    SCHUR_MAX_IT    = 500;
constexpr double SCHUR_TOL       = 1.0e-6;
constexpr int    VELOCITY_MAX_IT = 200;
constexpr double VELOCITY_TOL    = 1.0e-8;
constexpr int    AMG_V_CYCLES    = 2;

// ── Output ───────────────────────────────────────────────────────────────────
constexpr int         OUTPUT_INTERVAL = 5;
constexpr const char* OUTPUT_FILE     = "results/dg_channel_flow";

// ── Drag/lift normalization ──────────────────────────────────────────────────
constexpr double DRAG_LIFT_NORM = 2.0 / (RHO * U_MEAN * U_MEAN * CYL_DIAMETER);

} // namespace Params
