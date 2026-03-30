cmake_minimum_required(VERSION 3.16)
project(ChannelFlow LANGUAGES CXX)

# ── C++ standard ──────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)   # enables IDE / clangd support

# ── Locate libmesh-config ─────────────────────────────────────────────────────
#
# libmesh-config is a shell script installed with libMesh that reports all
# compiler and linker flags needed to build against it, including PETSc and
# any other dependencies libMesh was built with.
#
# Usage on the cluster (set LIBMESH_DIR before running cmake):
#   module load gcc openmpi petsc libmesh
#   cmake -DLIBMESH_DIR=$LIBMESH_DIR -B build .
#
find_program(LIBMESH_CONFIG
    NAMES libmesh-config
    HINTS "$ENV{LIBMESH_DIR}/bin" "${LIBMESH_DIR}/bin"
    DOC "libMesh configuration script (provides compile/link flags)")

if(NOT LIBMESH_CONFIG)
    message(FATAL_ERROR
        "libmesh-config not found.\n"
        "Set LIBMESH_DIR to the libMesh install prefix, e.g.:\n"
        "  cmake -DLIBMESH_DIR=/path/to/libmesh <srcdir>\n"
        "or load the cluster module before invoking cmake.")
endif()
message(STATUS "libmesh-config : ${LIBMESH_CONFIG}")

# ── Query libMesh build flags ─────────────────────────────────────────────────
#
# Helper: run libmesh-config with one flag, store result, abort on error.
#
function(libmesh_config_query flag outvar)
    execute_process(
        COMMAND "${LIBMESH_CONFIG}" "${flag}"
        OUTPUT_VARIABLE _val
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _rc)
    if(_rc)
        message(FATAL_ERROR
            "libmesh-config ${flag} failed (exit ${_rc}).\n"
            "Check that LIBMESH_DIR points to a valid libMesh installation.")
    endif()
    set(${outvar} "${_val}" PARENT_SCOPE)
endfunction()

libmesh_config_query(--cppflags LM_CPPFLAGS)   # preprocessor defines + -I paths
libmesh_config_query(--cxxflags LM_CXXFLAGS)   # C++ compiler flags (-std, -O2, …)
libmesh_config_query(--include  LM_INCLUDE)    # libMesh include directory (-I …)
libmesh_config_query(--libs     LM_LIBS)       # all link libraries and -L/-l flags

message(STATUS "libMesh cppflags: ${LM_CPPFLAGS}")
message(STATUS "libMesh cxxflags: ${LM_CXXFLAGS}")
message(STATUS "libMesh include : ${LM_INCLUDE}")

# ── Apply compile flags globally ──────────────────────────────────────────────
#
# All targets in this project are built against libMesh, so it is cleaner and
# more robust to apply the flags globally rather than per-target.
#
separate_arguments(_lm_cpp UNIX_COMMAND "${LM_CPPFLAGS}")
separate_arguments(_lm_cxx UNIX_COMMAND "${LM_CXXFLAGS}")
add_compile_options(${_lm_cpp} ${_lm_cxx})

# ── Apply include directories globally ───────────────────────────────────────
#
# Strip the leading -I from each flag to get a plain path for include_directories.
# Combine flags from both --include and --cppflags (the latter may also carry -I).
#
foreach(_raw_str "${LM_INCLUDE}" "${LM_CPPFLAGS}")
    string(REGEX MATCHALL "-I[^ ]+" _iflags "${_raw_str}")
    foreach(_f ${_iflags})
        string(SUBSTRING "${_f}" 2 -1 _path)
        include_directories("${_path}")
    endforeach()
endforeach()

# ── Store link flags for executable targets ───────────────────────────────────
#
# Using CMAKE_EXE_LINKER_FLAGS (a raw string, not a list) is the most portable
# approach for complex HPC link lines that include -Wl,… and RPATH flags.
#
string(APPEND CMAKE_EXE_LINKER_FLAGS " ${LM_LIBS}")

# ── Subdirectories ────────────────────────────────────────────────────────────
add_subdirectory(src)
enable_testing()
add_subdirectory(tests)
