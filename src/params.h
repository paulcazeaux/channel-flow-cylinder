#pragma once
/**
 * @file params.h
 * @brief Centralized numerical and physical parameters for the channel-flow solver.
 *
 * Every magic number in the C++ source must reference a constant defined here.
 * To change the Reynolds number, update U_MAX only; all derived quantities
 * (Re, U_MEAN, drag normalization) update automatically.
 *
 * Benchmark reference: Schafer & Turek (1996), DFG 2D-1.
 */

namespace Params {

// ── Physical domain (Schafer-Turek 2D-1) ─────────────────────────────────────
constexpr double CHANNEL_LENGTH = 2.2;   ///< Channel length [m]
constexpr double CHANNEL_HEIGHT = 0.41;  ///< Channel height H [m]
constexpr double CYL_X          = 0.2;   ///< Cylinder centre x [m]
constexpr double CYL_Y          = 0.2;   ///< Cylinder centre y [m]
constexpr double CYL_RADIUS     = 0.05;  ///< Cylinder radius r [m]
constexpr double CYL_DIAMETER   = 2.0 * CYL_RADIUS; ///< D = 2r [m]

// ── Fluid properties ──────────────────────────────────────────────────────────
constexpr double NU  = 1.0e-3;  ///< Kinematic viscosity ν [m²/s]
constexpr double RHO = 1.0;     ///< Fluid density ρ [kg/m³]

// ── Inlet velocity (parabolic profile) ───────────────────────────────────────
///   u(0, y) = 4 · U_MAX · y · (H − y) / H²,   v = 0
///
/// U_MAX = 0.3  m/s  →  U_MEAN = 0.2  m/s  →  Re = 20  (Schafer-Turek ref.)
/// U_MAX = 0.075 m/s →  U_MEAN = 0.05 m/s  →  Re =  5  (primary target)
constexpr double U_MAX  = 0.3;                    ///< Peak inlet velocity [m/s]
constexpr double U_MEAN = (2.0 / 3.0) * U_MAX;   ///< Mean inlet velocity [m/s]

/// Reynolds number: Re = U_MEAN · D / ν
constexpr double RE = U_MEAN * CYL_DIAMETER / NU;

// ── Boundary IDs — must match Physical Curve tags in channel.geo ─────────────
constexpr int BID_INLET    = 1;  ///< Inlet   (x = 0)
constexpr int BID_OUTLET   = 2;  ///< Outlet  (x = CHANNEL_LENGTH)
constexpr int BID_WALLS    = 3;  ///< Top and bottom walls
constexpr int BID_CYLINDER = 4;  ///< Cylinder surface

// ── FEM discretisation ────────────────────────────────────────────────────────
constexpr int VELOCITY_ORDER = 2;  ///< Taylor-Hood: P2 velocity
constexpr int PRESSURE_ORDER = 1;  ///< Taylor-Hood: P1 pressure

// ── Nonlinear solver: PETSc SNES (Newton + backtracking line search) ──────────
constexpr int    SNES_MAX_IT = 50;     ///< Maximum Newton iterations
constexpr double SNES_ATOL   = 1.0e-8; ///< Absolute nonlinear residual tolerance
constexpr double SNES_RTOL   = 1.0e-10;///< Relative nonlinear residual tolerance

// ── Linear solver: FGMRES + BoomerAMG (per Newton step) ──────────────────────
constexpr int    KSP_MAX_IT    = 500;    ///< Maximum Krylov iterations
constexpr double KSP_RTOL      = 1.0e-10;///< Relative linear residual tolerance
constexpr int    GMRES_RESTART = 100;    ///< FGMRES restart parameter

// ── Drag/lift normalisation ───────────────────────────────────────────────────
///   C_D = 2 F_D / (ρ U_MEAN² D),   C_L = 2 F_L / (ρ U_MEAN² D)
constexpr double DRAG_LIFT_NORM = 2.0 / (RHO * U_MEAN * U_MEAN * CYL_DIAMETER);

// ── Output ────────────────────────────────────────────────────────────────────
constexpr const char* OUTPUT_FILE = "results/channel_flow.e"; ///< ExodusII output path

} // namespace Params
