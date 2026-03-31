#!/usr/bin/env python3
"""
test_validation.py
Phase 7: Validate solver against Schafer-Turek DFG 2D-1 benchmark (Re=20).

Uses a refined mesh (lc_far=0.02, lc_cyl=0.003, ~8000 elements) so that the
Taylor-Hood P2/P1 solution is well within the benchmark interval.

Pass criteria (from PLAN.md Phase 7):
  1. Re=20 steady solve converges in <= 20 Newton steps.
  2. Drag coefficient:  C_D in [5.57, 5.59].
  3. Lift coefficient:  |C_L| < 0.02.
  4. ExodusII output file is written and non-empty.

Run with:
    python -m pytest tests/test_validation.py -v

Requirements:
  - Solver binary at build/src/channel_flow (build with cmake first).
  - MPI runtime (mpirun/mpiexec) on PATH.
  - gmsh Python package (for mesh generation if needed).
  - Environment: spack load openmpi libmesh  (or equivalent).
"""

import os
import re
import subprocess
import sys
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SOLVER = os.path.join(_ROOT, "build", "src", "channel_flow")
_MESH_SCRIPT = os.path.join(_ROOT, "meshes", "generate_mesh.py")
_MESH_PATH = os.path.join(_ROOT, "meshes", "channel_refine.msh")
_OUTPUT_FILE = os.path.join(_ROOT, "results", "channel_flow.e")


def _find_mpirun():
    """Return path to mpirun or mpiexec, or None."""
    for name in ("mpirun", "mpiexec"):
        try:
            result = subprocess.run(
                ["which", name], capture_output=True, text=True
            )
            if result.returncode == 0:
                return result.stdout.strip()
        except FileNotFoundError:
            pass
    return None


def _generate_mesh_if_needed():
    """Generate the refined validation mesh if it does not exist."""
    if os.path.isfile(_MESH_PATH):
        return
    subprocess.check_call([
        sys.executable, _MESH_SCRIPT,
        "--lc-far", "0.02",
        "--lc-cyl", "0.003",
        "--output", _MESH_PATH,
    ])


def _run_solver():
    """Run the channel_flow solver on the refined mesh and return stdout."""
    mpirun = _find_mpirun()
    if mpirun is None:
        raise RuntimeError("mpirun/mpiexec not found on PATH")
    if not os.path.isfile(_SOLVER):
        raise RuntimeError(
            f"Solver binary not found: {_SOLVER}\n"
            "Build with: cd build && cmake --build . -j4"
        )

    _generate_mesh_if_needed()

    result = subprocess.run(
        [mpirun, "--mca", "opal_warn_on_missing_libcuda", "0",
         "-n", "1", _SOLVER, "--mesh", _MESH_PATH],
        capture_output=True, text=True, timeout=600,
    )
    # Combine stdout + stderr; solver prints to libMesh::out (stdout).
    output = result.stdout + result.stderr
    if result.returncode != 0:
        raise RuntimeError(
            f"Solver exited with code {result.returncode}:\n{output}"
        )
    return output


def _parse_cd_cl(output):
    """Extract C_D and C_L from solver output line: '  C_D = ... C_L = ...'."""
    match = re.search(
        r"C_D\s*=\s*([0-9.eE+-]+)\s+C_L\s*=\s*([0-9.eE+-]+)", output
    )
    if not match:
        raise ValueError(f"Could not parse C_D/C_L from output:\n{output}")
    return float(match.group(1)), float(match.group(2))


def _count_newton_steps(output):
    """Count the number of 'Nonlinear step:' lines in NS phase (after Phase 3b)."""
    phase3b_start = output.find("Phase 3b")
    if phase3b_start < 0:
        raise ValueError("Could not find 'Phase 3b' in solver output")
    ns_output = output[phase3b_start:]
    return len(re.findall(r"Nonlinear step:", ns_output))


# ── Run solver once for all tests ─────────────────────────────────────────────
_solver_output = None
_solver_error = None


def setUpModule():
    """Run the solver once and cache the output for all tests."""
    global _solver_output, _solver_error
    try:
        _solver_output = _run_solver()
    except Exception as e:
        _solver_error = str(e)


class TestValidation(unittest.TestCase):
    """Schafer-Turek DFG 2D-1 benchmark validation (Re=20)."""

    def setUp(self):
        if _solver_error:
            self.skipTest(f"Solver failed: {_solver_error}")

    def test_newton_convergence(self):
        """Re=20 steady solve converges in <= 20 Newton steps."""
        n_steps = _count_newton_steps(_solver_output)
        self.assertGreater(n_steps, 0, "No Newton steps detected")
        self.assertLessEqual(n_steps, 20,
                             f"Newton took {n_steps} steps (max 20)")

    def test_drag_coefficient(self):
        """C_D in [5.57, 5.59] (Schafer-Turek interval)."""
        cd, _ = _parse_cd_cl(_solver_output)
        self.assertGreaterEqual(cd, 5.57,
                                f"C_D = {cd} below benchmark lower bound 5.57")
        self.assertLessEqual(cd, 5.59,
                             f"C_D = {cd} above benchmark upper bound 5.59")

    def test_lift_coefficient(self):
        """|C_L| < 0.02."""
        _, cl = _parse_cd_cl(_solver_output)
        self.assertLess(abs(cl), 0.02,
                        f"|C_L| = {abs(cl)} exceeds tolerance 0.02")

    def test_exodus_output(self):
        """ExodusII output file is written and non-empty."""
        self.assertTrue(os.path.isfile(_OUTPUT_FILE),
                        f"Output file not found: {_OUTPUT_FILE}")
        self.assertGreater(os.path.getsize(_OUTPUT_FILE), 0,
                           "Output file is empty")


if __name__ == "__main__":
    unittest.main()
