#!/usr/bin/env python3
"""
test_driver.py
Tests for Python driver CLI parsing, SLURM generation, mesh invocation.

Pass criteria:
  1. CLI argument parsing works — all args parsed without error.
  2. SLURM script has correct structure — #SBATCH directives present.
  3. Mesh generation is invoked correctly — generate_mesh.py called with right args.
  4. Dry-run mode does not submit or execute — no side effects.
  5. Solver options (--re, --dt, --t-final, --output-interval) are parsed.

Run with:
    python -m pytest tests/test_driver.py -v
or:
    python tests/test_driver.py
"""

import os
import sys
import unittest

# Add scripts/ to path so we can import run_simulation
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

import run_simulation


class TestCLIParsing(unittest.TestCase):
    """CLI argument parsing works — all args parsed without error."""

    def test_defaults(self):
        """Default arguments are parsed correctly."""
        args = run_simulation.parse_args([])
        self.assertIsNone(args.mesh_file)
        self.assertEqual(args.lc_far, 0.05)
        self.assertEqual(args.lc_cyl, 0.01)
        self.assertEqual(args.re, 10.0)
        self.assertEqual(args.dt, 0.025)
        self.assertEqual(args.t_final, 8.0)
        self.assertEqual(args.output_interval, 10)
        self.assertEqual(args.nodes, 1)
        self.assertEqual(args.ntasks, 1)
        self.assertEqual(args.time, "01:00:00")
        self.assertEqual(args.partition, "batch")
        self.assertEqual(args.job_name, "channel_flow")
        self.assertEqual(args.output_dir, "results")
        self.assertFalse(args.dry_run)

    def test_custom_args(self):
        """Custom arguments override defaults."""
        args = run_simulation.parse_args([
            "--lc-far", "0.1",
            "--lc-cyl", "0.03",
            "--re", "20",
            "--dt", "0.01",
            "--t-final", "5.0",
            "--output-interval", "5",
            "--nodes", "2",
            "--ntasks", "4",
            "--time", "02:00:00",
            "--partition", "gpu",
            "--job-name", "my_job",
            "--output-dir", "/tmp/out",
            "--dry-run",
        ])
        self.assertEqual(args.lc_far, 0.1)
        self.assertEqual(args.lc_cyl, 0.03)
        self.assertEqual(args.re, 20.0)
        self.assertEqual(args.dt, 0.01)
        self.assertEqual(args.t_final, 5.0)
        self.assertEqual(args.output_interval, 5)
        self.assertEqual(args.nodes, 2)
        self.assertEqual(args.ntasks, 4)
        self.assertEqual(args.time, "02:00:00")
        self.assertEqual(args.partition, "gpu")
        self.assertEqual(args.job_name, "my_job")
        self.assertEqual(args.output_dir, "/tmp/out")
        self.assertTrue(args.dry_run)

    def test_mesh_file_arg(self):
        """--mesh-file is parsed correctly."""
        args = run_simulation.parse_args(["--mesh-file", "/path/to/mesh.msh"])
        self.assertEqual(args.mesh_file, "/path/to/mesh.msh")

    def test_solver_args(self):
        """Solver-specific args are parsed correctly."""
        args = run_simulation.parse_args([
            "--re", "5",
            "--dt", "0.05",
            "--t-final", "12.0",
            "--output-interval", "20",
        ])
        self.assertEqual(args.re, 5.0)
        self.assertEqual(args.dt, 0.05)
        self.assertEqual(args.t_final, 12.0)
        self.assertEqual(args.output_interval, 20)


class TestSLURMScript(unittest.TestCase):
    """SLURM script is generated with correct structure."""

    def setUp(self):
        self.args = run_simulation.parse_args([
            "--nodes", "2",
            "--ntasks", "8",
            "--time", "04:00:00",
            "--partition", "debug",
            "--job-name", "test_job",
            "--output-dir", "test_results",
        ])
        self.script = run_simulation.generate_slurm_script(
            self.args, "/path/to/mesh.msh"
        )

    def test_shebang(self):
        """Script starts with bash shebang."""
        self.assertTrue(self.script.startswith("#!/bin/bash"))

    def test_sbatch_directives(self):
        """All #SBATCH directives are present."""
        self.assertIn("#SBATCH --job-name=test_job", self.script)
        self.assertIn("#SBATCH --nodes=2", self.script)
        self.assertIn("#SBATCH --ntasks=8", self.script)
        self.assertIn("#SBATCH --time=04:00:00", self.script)
        self.assertIn("#SBATCH --partition=debug", self.script)
        self.assertIn("#SBATCH --output=test_job_%j.out", self.script)
        self.assertIn("#SBATCH --error=test_job_%j.err", self.script)

    def test_module_loads(self):
        """Module load commands are present."""
        self.assertIn("module load", self.script)

    def test_mpirun(self):
        """mpirun command uses correct ntasks and mesh path."""
        self.assertIn("mpirun -n 8", self.script)
        self.assertIn("--mesh /path/to/mesh.msh", self.script)

    def test_mkdir_output(self):
        """Output directory is created."""
        self.assertIn("mkdir -p test_results", self.script)


class TestMeshGeneration(unittest.TestCase):
    """Mesh generation is invoked correctly."""

    def test_existing_mesh_skips_generation(self):
        """When --mesh-file is provided, generate_mesh returns it directly."""
        args = run_simulation.parse_args([
            "--mesh-file", "/existing/mesh.msh",
            "--dry-run",
        ])
        result = run_simulation.generate_mesh(args)
        self.assertEqual(result, "/existing/mesh.msh")


class TestDryRun(unittest.TestCase):
    """Dry-run mode does not submit or execute."""

    def test_dry_run_no_side_effects(self):
        """Dry-run completes without creating files or submitting jobs."""
        _, mesh_path, script = run_simulation.main([
            "--mesh-file", "/tmp/fake.msh",
            "--dry-run",
            "--job-name", "dryrun_test",
        ])
        self.assertEqual(mesh_path, "/tmp/fake.msh")
        self.assertIn("#SBATCH", script)

        # Verify no SLURM script file was written
        slurm_path = os.path.join(_ROOT, "dryrun_test.slurm")
        self.assertFalse(os.path.exists(slurm_path))


if __name__ == "__main__":
    unittest.main()
