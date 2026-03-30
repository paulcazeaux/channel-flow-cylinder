# tests/CMakeLists.txt
# CTest registration for C++ unit tests and mesh generation utilities.
#
# Run all tests:   ctest --output-on-failure  (from the build directory)
# Run Python tests: python -m pytest <srcdir>/tests/ -v

find_package(Python3 COMPONENTS Interpreter REQUIRED)

# ── Coarse test mesh ──────────────────────────────────────────────────────────
#
# A coarse mesh (lc_far=0.1, lc_cyl=0.03) is generated on demand and used by
# all C++ tests so they complete in seconds on a workstation.
# The file is git-ignored; regenerate with: cmake --build . --target generate_coarse_mesh
#
set(COARSE_MESH "${CMAKE_SOURCE_DIR}/meshes/channel_coarse.msh")

add_custom_target(generate_coarse_mesh
    COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_SOURCE_DIR}/meshes/generate_mesh.py"
            --lc-far 0.1 --lc-cyl 0.03
            --output "${COARSE_MESH}"
    COMMENT "Generating coarse test mesh (lc_far=0.1, lc_cyl=0.03)"
    BYPRODUCTS "${COARSE_MESH}")

# ── Helper macro: build and register a C++ unit test ─────────────────────────
#
# Each test executable is linked against channel_flow_lib and receives the
# coarse mesh path via the TEST_MESH environment variable.
#
macro(add_solver_test name source)
    add_executable(${name} ${source})
    target_include_directories(${name} PRIVATE "${CMAKE_SOURCE_DIR}/src")
    target_link_libraries(${name} PRIVATE channel_flow_lib)
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES
        ENVIRONMENT "TEST_MESH=${COARSE_MESH}"
        TIMEOUT 120)
endmacro()

# ── Phase 2: FEM space initialisation ────────────────────────────────────────
add_solver_test(test_spaces test_spaces.cpp)

# ── Phase 3: Stokes linear solve ─────────────────────────────────────────────
# add_solver_test(test_stokes test_stokes.cpp)   # uncomment after Phase 3

# ── Phase 4: ExodusII output + drag/lift ─────────────────────────────────────
# add_solver_test(test_output test_output.cpp)   # uncomment after Phase 4
