# src/CMakeLists.txt
# Builds the channel-flow solver library and main executable.
#
# Target layout:
#   channel_flow_lib  STATIC  — solver components (no main); linked by tests
#   channel_flow      EXE     — main entry point, links channel_flow_lib

# ── Solver source files ───────────────────────────────────────────────────────
set(SOLVER_SOURCES
    channel_flow_system.cpp    # DG variable setup, context building, fieldsplit
    channel_flow_assembly.cpp  # volume (element-interior) assembly
    dg_face_assembly.cpp       # DG face integrals (SIPG, upwind, pressure)
    drag_lift.cpp              # stress-traction boundary integration
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
