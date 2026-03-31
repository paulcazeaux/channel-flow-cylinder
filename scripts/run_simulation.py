#!/usr/bin/env python3
"""
run_simulation.py
Python driver for the channel-flow-cylinder solver.

Orchestrates mesh generation, SLURM job script creation, and (optionally)
job submission.  Designed to run on the Oscar HPC cluster at Brown.

Usage:
    python scripts/run_simulation.py [OPTIONS]

    --mesh-file PATH     Path to an existing .msh mesh (skip generation)
    --lc-far FLOAT       Global mesh size [m]      (default: 0.05)
    --lc-cyl FLOAT       Cylinder mesh size [m]    (default: 0.01)
    --nodes INT           SLURM nodes               (default: 1)
    --ntasks INT          SLURM MPI tasks            (default: 1)
    --time STR            SLURM wall-clock limit     (default: 01:00:00)
    --partition STR       SLURM partition            (default: batch)
    --job-name STR        SLURM job name             (default: channel_flow)
    --dry-run             Print actions without executing
    --output-dir PATH     Directory for results      (default: results)
"""

import argparse
import os
import subprocess
import sys
import textwrap


# ── Project layout ────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_SCRIPT_DIR)
_MESH_SCRIPT = os.path.join(_ROOT, "meshes", "generate_mesh.py")
_DEFAULT_MESH = os.path.join(_ROOT, "meshes", "channel.msh")
_BUILD_DIR = os.path.join(_ROOT, "build")
_SOLVER_BIN = os.path.join(_BUILD_DIR, "src", "channel_flow")


def parse_args(argv=None):
    """Parse command-line arguments.

    Args:
        argv: Argument list (defaults to sys.argv[1:]).

    Returns:
        argparse.Namespace with all options.
    """
    p = argparse.ArgumentParser(
        description="Driver for 2D channel-flow-cylinder solver.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Mesh options
    mesh = p.add_argument_group("mesh")
    mesh.add_argument("--mesh-file", type=str, default=None,
                      help="Path to existing .msh file (skip generation)")
    mesh.add_argument("--lc-far", type=float, default=0.05,
                      help="Global element size [m]")
    mesh.add_argument("--lc-cyl", type=float, default=0.01,
                      help="Cylinder element size [m]")

    # SLURM options
    slurm = p.add_argument_group("SLURM")
    slurm.add_argument("--nodes", type=int, default=1,
                       help="Number of nodes")
    slurm.add_argument("--ntasks", type=int, default=1,
                       help="Number of MPI tasks")
    slurm.add_argument("--time", type=str, default="01:00:00",
                       help="Wall-clock time limit")
    slurm.add_argument("--partition", type=str, default="batch",
                       help="SLURM partition")
    slurm.add_argument("--job-name", type=str, default="channel_flow",
                       help="SLURM job name")

    # Output / control
    p.add_argument("--output-dir", type=str, default="results",
                   help="Directory for solver output")
    p.add_argument("--dry-run", action="store_true",
                   help="Print actions without executing")

    return p.parse_args(argv)


def generate_mesh(args):
    """Invoke meshes/generate_mesh.py to create the mesh.

    Args:
        args: Parsed arguments (uses lc_far, lc_cyl, mesh_file, dry_run).

    Returns:
        Path to the mesh file.
    """
    if args.mesh_file:
        return args.mesh_file

    mesh_path = _DEFAULT_MESH
    cmd = [
        sys.executable, _MESH_SCRIPT,
        "--lc-far", str(args.lc_far),
        "--lc-cyl", str(args.lc_cyl),
        "--output", mesh_path,
    ]

    print(f"[driver] Generating mesh: {' '.join(cmd)}")
    if not args.dry_run:
        subprocess.check_call(cmd)

    return mesh_path


def generate_slurm_script(args, mesh_path):
    """Build a SLURM batch script as a string.

    Args:
        args: Parsed arguments (SLURM options).
        mesh_path: Path to the mesh file.

    Returns:
        SLURM script content as a string.
    """
    solver = _SOLVER_BIN
    return textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name={args.job_name}
        #SBATCH --nodes={args.nodes}
        #SBATCH --ntasks={args.ntasks}
        #SBATCH --time={args.time}
        #SBATCH --partition={args.partition}
        #SBATCH --output={args.job_name}_%j.out
        #SBATCH --error={args.job_name}_%j.err

        module load gcc openmpi petsc libmesh

        mkdir -p {args.output_dir}
        mpirun -n {args.ntasks} {solver} --mesh {mesh_path}
    """)


def main(argv=None):
    """Entry point: parse args, generate mesh, write SLURM script, optionally submit."""
    args = parse_args(argv)

    # Step 1: mesh
    mesh_path = generate_mesh(args)

    # Step 2: SLURM script
    script = generate_slurm_script(args, mesh_path)
    script_path = os.path.join(_ROOT, f"{args.job_name}.slurm")

    print(f"[driver] Writing SLURM script: {script_path}")
    if not args.dry_run:
        with open(script_path, "w") as f:
            f.write(script)

    # Step 3: submit (unless dry-run)
    if args.dry_run:
        print("[driver] Dry-run mode — no submission.")
        print("--- SLURM script ---")
        print(script)
    else:
        print(f"[driver] Submitting: sbatch {script_path}")
        subprocess.check_call(["sbatch", script_path])

    return args, mesh_path, script


if __name__ == "__main__":
    main()
