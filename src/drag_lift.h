#pragma once
/**
 * @file drag_lift.h
 * @brief Drag and lift force computation by boundary integration over BID_CYLINDER.
 *
 * The stress traction is integrated over the cylinder boundary using the
 * Cauchy stress tensor for a Newtonian fluid:
 *
 *   σ = -p I + ν (∇u + ∇uᵀ)
 *
 * The drag and lift forces are the x- and y-components of ∫_Γ σ·n dS.
 * Normalize with Params::DRAG_LIFT_NORM to obtain C_D and C_L.
 */

#include "libmesh/mesh_base.h"

#include <utility>

class ChannelFlowSystem;

/**
 * @brief Integrate the viscous stress tensor over the cylinder surface.
 *
 * Iterates over all sides of active local elements tagged with
 * Params::BID_CYLINDER, performs side FE quadrature (P2 velocity, P1 pressure),
 * and accumulates the x/y traction components.  Reduces the local sums across
 * all MPI ranks before returning.
 *
 * Must be called after sys.solve() so that sys.current_local_solution is current.
 *
 * @param sys  Solved ChannelFlowSystem.
 * @param mesh The mesh (queried for BID_CYLINDER boundary tags).
 * @return {F_D, F_L} — raw forces in Newtons (before normalization).
 */
std::pair<double, double> compute_drag_lift(const ChannelFlowSystem& sys,
                                            const libMesh::MeshBase& mesh);
