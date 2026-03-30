# src/CMakeLists.txt
# Builds the channel-flow solver library and main executable.
#
# Target layout:
#   channel_flow_lib  STATIC  — solver components (no main); linked by tests
#   channel_flow      EXE     — main entry point, links channel_flow_lib
#
# Phase 5: SOLVER_SOURCES is empty; channel_flow_lib is not created yet.
# Phase 2: add channel_flow_system.cpp to SOLVER_SOURCES and uncomment the lib.

# ── Solver source files (populated in Phase 2) ───────────────────────────────
set(SOLVER_SOURCES
    # channel_flow_system.cpp   # Phase 2: FEM assembly and BCs
)

# ── channel_flow_lib: reusable solver components ─────────────────────────────
#
# This library is the link target for the C++ unit tests (test_spaces,
# test_stokes, test_output).  Created only when there are sources to compile.
#
if(SOLVER_SOURCES)
    add_library(channel_flow_lib STATIC ${SOLVER_SOURCES})
    target_include_directories(channel_flow_lib
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endif()

# ── channel_flow: main executable ────────────────────────────────────────────
add_executable(channel_flow main.cpp)
target_include_directories(channel_flow
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
if(SOLVER_SOURCES)
    target_link_libraries(channel_flow PRIVATE channel_flow_lib)
endif()

install(TARGETS channel_flow DESTINATION bin)
