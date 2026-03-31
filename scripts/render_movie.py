#!/usr/bin/env pvpython
"""
render_movie.py
Render velocity magnitude frames from an ExodusII time series using ParaView.

Produces PNG frames in results/frames/, then optionally assembles into an
MP4 movie using ffmpeg.

Usage:
    pvpython scripts/render_movie.py [--input FILE] [--output-dir DIR]
                                     [--width W] [--height H] [--no-ffmpeg]

Requires: ParaView (pvpython), optionally ffmpeg for movie assembly.
"""

import argparse
import os
import sys

def parse_args():
    p = argparse.ArgumentParser(
        description="Render velocity magnitude movie from ExodusII file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", type=str, default="results/channel_flow.e",
                   help="Path to ExodusII file")
    p.add_argument("--output-dir", type=str, default="results/frames",
                   help="Directory for PNG frames")
    p.add_argument("--movie", type=str, default="results/channel_flow.mp4",
                   help="Output movie path")
    p.add_argument("--width", type=int, default=1920,
                   help="Frame width in pixels")
    p.add_argument("--height", type=int, default=480,
                   help="Frame height (channel is 5:1 aspect)")
    p.add_argument("--no-ffmpeg", action="store_true",
                   help="Skip ffmpeg movie assembly")
    p.add_argument("--colormap", type=str, default="Cool to Warm",
                   help="ParaView colormap name")
    p.add_argument("--field", type=str, default="velocity_magnitude",
                   choices=["velocity_magnitude", "u", "v", "p"],
                   help="Field to visualise")
    return p.parse_args()


def render(args):
    from paraview.simple import (
        ExodusIIReader, GetActiveViewOrCreate, Show, ColorBy,
        GetColorTransferFunction, GetOpacityTransferFunction,
        GetScalarBar, HideScalarBarIfNotNeeded, SaveScreenshot,
        GetAnimationScene, Calculator, GetDisplayProperties,
    )

    os.makedirs(args.output_dir, exist_ok=True)

    # ── Load data ──────────────────────────────────────────────────────────
    reader = ExodusIIReader(FileName=args.input)
    reader.UpdatePipeline()

    # Select point variables
    reader.PointVariables = ['u', 'v', 'p']
    reader.UpdatePipeline()

    # ── Compute velocity magnitude ─────────────────────────────────────────
    if args.field == "velocity_magnitude":
        calc = Calculator(Input=reader)
        calc.Function = 'sqrt(u^2 + v^2)'
        calc.ResultArrayName = 'vel_mag'
        calc.UpdatePipeline()
        source = calc
        color_field = ('POINTS', 'vel_mag')
    else:
        source = reader
        color_field = ('POINTS', args.field)

    # ── Setup view ─────────────────────────────────────────────────────────
    view = GetActiveViewOrCreate('RenderView')
    view.ViewSize = [args.width, args.height]
    view.Background = [1.0, 1.0, 1.0]  # white background

    display = Show(source, view)
    ColorBy(display, color_field)

    # Camera: top-down view of the 2D domain
    view.CameraPosition = [1.1, 0.205, 3.0]
    view.CameraFocalPoint = [1.1, 0.205, 0.0]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.CameraParallelProjection = 1
    view.CameraParallelScale = 0.25

    # Colormap
    field_name = 'vel_mag' if args.field == 'velocity_magnitude' else args.field
    lut = GetColorTransferFunction(field_name)
    lut.ApplyPreset(args.colormap, True)

    # Set range based on field
    if args.field == "velocity_magnitude":
        lut.RescaleTransferFunction(0.0, 2.0)
    elif args.field == "p":
        lut.RescaleTransferFunction(-0.5, 0.5)
    else:
        lut.RescaleTransferFunction(-0.5, 2.0)

    # Show color bar
    display.SetScalarBarVisibility(view, True)
    bar = GetScalarBar(lut, view)
    bar.Title = field_name
    bar.ComponentTitle = ''

    # ── Render frames ──────────────────────────────────────────────────────
    scene = GetAnimationScene()
    timesteps = reader.TimestepValues

    print(f"Rendering {len(timesteps)} frames to {args.output_dir}/")
    for i, t in enumerate(timesteps):
        scene.AnimationTime = t
        view.Update()

        frame_path = os.path.join(args.output_dir, f"frame_{i:04d}.png")
        SaveScreenshot(frame_path, view,
                       ImageResolution=[args.width, args.height])

        if (i + 1) % 10 == 0 or i == len(timesteps) - 1:
            print(f"  Frame {i+1}/{len(timesteps)}  t={t:.3f}")

    print(f"Done. {len(timesteps)} frames in {args.output_dir}/")
    return len(timesteps)


def assemble_movie(args, n_frames):
    """Assemble PNG frames into MP4 using ffmpeg."""
    import subprocess

    cmd = [
        "ffmpeg", "-y",
        "-framerate", "15",
        "-i", os.path.join(args.output_dir, "frame_%04d.png"),
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        "-crf", "18",
        args.movie,
    ]

    print(f"Assembling movie: {' '.join(cmd)}")
    try:
        subprocess.check_call(cmd)
        print(f"Movie written to {args.movie}")
    except FileNotFoundError:
        print("ffmpeg not found. Install it or use --no-ffmpeg.")
        print(f"Frames are in {args.output_dir}/ — assemble manually with:")
        print(f"  ffmpeg -framerate 15 -i {args.output_dir}/frame_%04d.png "
              f"-c:v libx264 -pix_fmt yuv420p {args.movie}")
    except subprocess.CalledProcessError as e:
        print(f"ffmpeg failed: {e}")


def main():
    args = parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: input file not found: {args.input}")
        sys.exit(1)

    n_frames = render(args)

    if not args.no_ffmpeg:
        assemble_movie(args, n_frames)


if __name__ == "__main__":
    main()
