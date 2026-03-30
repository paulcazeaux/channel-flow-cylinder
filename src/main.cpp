/**
 * @file main.cpp
 * @brief Entry point for the channel-flow solver.
 *
 * Phase 5: minimal stub that initialises libMesh and prints the problem
 * configuration derived from params.h.  This validates that the build system
 * correctly finds libMesh headers and links the library.
 *
 * Phase 2 will replace this stub with the full Newton solver loop.
 */

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"

#include "params.h"

int main(int argc, char** argv)
{
    // Initialise libMesh (and MPI if available).
    // LibMeshInit must be the first object constructed and last destroyed.
    libMesh::LibMeshInit init(argc, argv);

    libMesh::out
        << "ChannelFlow: 2D incompressible Navier-Stokes past a cylinder\n"
        << "  Domain  : [0, " << Params::CHANNEL_LENGTH << "] x "
                              << "[0, " << Params::CHANNEL_HEIGHT << "] m\n"
        << "  Cylinder: centre (" << Params::CYL_X << ", " << Params::CYL_Y
                                  << ")  r = " << Params::CYL_RADIUS << " m\n"
        << "  U_MAX   : " << Params::U_MAX   << " m/s\n"
        << "  U_MEAN  : " << Params::U_MEAN  << " m/s\n"
        << "  Re      : " << Params::RE      << "\n"
        << "  Elements: P" << Params::VELOCITY_ORDER
                           << "/P" << Params::PRESSURE_ORDER
                           << " (Taylor-Hood)\n"
        << "  SNES tol: atol=" << Params::SNES_ATOL
                               << "  rtol=" << Params::SNES_RTOL << "\n"
        << "  Output  : " << Params::OUTPUT_FILE << "\n"
        << "\n"
        << "  (Phase 5 build stub — solver not yet implemented)\n";

    return 0;
}
