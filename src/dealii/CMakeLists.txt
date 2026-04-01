# src/dealii/CMakeLists.txt
# DG Navier-Stokes solver using deal.II.
#
# Build:
#   cmake -DDEAL_II_DIR=<path> -B build_dealii src/dealii
#   cmake --build build_dealii

cmake_minimum_required(VERSION 3.13)

find_package(deal.II 9.7 REQUIRED
    HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR})

deal_ii_initialize_cached_variables()

project(dg_channel_flow CXX)

# ── Solver library ───────────────────────────────────────────────────────────
add_library(dg_solver_lib STATIC
    dg_navier_stokes.cc
    dg_assembly.cc
    dg_schur_solver.cc
)
deal_ii_setup_target(dg_solver_lib)
target_include_directories(dg_solver_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# ── Main executable ──────────────────────────────────────────────────────────
add_executable(dg_channel_flow main.cc)
deal_ii_setup_target(dg_channel_flow)
target_link_libraries(dg_channel_flow dg_solver_lib)

install(TARGETS dg_channel_flow DESTINATION bin)
