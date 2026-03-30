# tests/CMakeLists.txt
# CTest registration for C++ unit tests.
# Each test is a small standalone executable linked against libMesh + PETSc.
# Add tests here as phases are implemented; run with: ctest --output-on-failure

# ── Helper macro: build a test executable and register it with CTest ──────────
macro(add_solver_test name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE channel_flow_lib)
    add_test(NAME ${name} COMMAND ${name})
    # Use a coarse mesh so tests complete quickly on a workstation
    set_tests_properties(${name} PROPERTIES
        ENVIRONMENT "TEST_MESH=${CMAKE_SOURCE_DIR}/../meshes/channel_coarse.msh"
        TIMEOUT 120
    )
endmacro()

# ── Phase 2: FEM space initialisation ────────────────────────────────────────
# add_solver_test(test_spaces test_spaces.cpp)   # uncomment after Phase 2

# ── Phase 3: Stokes linear solve ─────────────────────────────────────────────
# add_solver_test(test_stokes test_stokes.cpp)   # uncomment after Phase 3

# ── Phase 4: ExodusII output + drag/lift ─────────────────────────────────────
# add_solver_test(test_output test_output.cpp)   # uncomment after Phase 4
