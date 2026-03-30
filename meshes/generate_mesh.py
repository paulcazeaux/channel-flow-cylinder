#!/usr/bin/env python3
"""
generate_mesh.py
Generate the triangular mesh for 2D channel flow past a cylinder.

Builds the Schafer-Turek 2D-1 domain using the Gmsh Python API,
applies a smooth refinement field near the cylinder, and exports
the mesh in Gmsh format 2.2 ASCII, which libMesh GmshIO can read.

Physical boundary tags written to the mesh (match boundary_id in main.cpp):
    1 = inlet    (x = 0)
    2 = outlet   (x = 2.2 m)
    3 = walls    (y = 0 and y = 0.41 m)
    4 = cylinder surface

Usage:
    python generate_mesh.py [--lc-far F] [--lc-cyl C] [--output FILE]

    --lc-far  Global characteristic element size in metres (default: 0.05)
    --lc-cyl  Near-cylinder element size in metres         (default: 0.01)
    --output  Output .msh file path                        (default: channel.msh
              next to this script)
"""

import argparse
import os
import sys

try:
    import gmsh
except ImportError:
    sys.exit(
        "Error: gmsh Python package not found.\n"
        "Install with:  pip install gmsh\n"
        "Or load the cluster module that provides it."
    )


# ── Domain constants (Schafer-Turek 2D-1) ────────────────────────────────────
_H  = 0.41   # channel height [m]
_L  = 2.2    # channel length [m]
_XC = 0.2    # cylinder center x [m]
_YC = 0.2    # cylinder center y [m]
_R  = 0.05   # cylinder radius [m]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate Gmsh mesh for channel-cylinder domain.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--lc-far", type=float, default=0.05,
                   metavar="F", help="Global characteristic element size [m]")
    p.add_argument("--lc-cyl", type=float, default=0.01,
                   metavar="C", help="Near-cylinder element size [m]")
    p.add_argument("--output", type=str,
                   default=os.path.join(os.path.dirname(__file__), "channel.msh"),
                   metavar="FILE", help="Output mesh file path")
    return p.parse_args()


def build_geometry(lc_far: float, lc_cyl: float) -> tuple:
    """
    Define the channel-cylinder geometry using the Gmsh built-in kernel.

    Returns:
        (surf, arcs) — surface tag and list of 4 cylinder arc tags,
        needed to set up the refinement field and physical groups.
    """
    # ── Channel corners ───────────────────────────────────────────────────────
    p1 = gmsh.model.geo.addPoint(0,    0,    0, lc_far)  # bottom-left
    p2 = gmsh.model.geo.addPoint(_L,   0,    0, lc_far)  # bottom-right
    p3 = gmsh.model.geo.addPoint(_L,   _H,   0, lc_far)  # top-right
    p4 = gmsh.model.geo.addPoint(0,    _H,   0, lc_far)  # top-left

    # ── Channel edges ─────────────────────────────────────────────────────────
    l_bottom = gmsh.model.geo.addLine(p1, p2)   # bottom wall
    l_outlet = gmsh.model.geo.addLine(p2, p3)   # outlet
    l_top    = gmsh.model.geo.addLine(p3, p4)   # top wall
    l_inlet  = gmsh.model.geo.addLine(p4, p1)   # inlet

    # ── Cylinder: center point + 4 boundary points ────────────────────────────
    pc = gmsh.model.geo.addPoint(_XC,      _YC,      0, lc_cyl)  # center (arc pivot)
    pe = gmsh.model.geo.addPoint(_XC + _R, _YC,      0, lc_cyl)  # east
    pn = gmsh.model.geo.addPoint(_XC,      _YC + _R, 0, lc_cyl)  # north
    pw = gmsh.model.geo.addPoint(_XC - _R, _YC,      0, lc_cyl)  # west
    ps = gmsh.model.geo.addPoint(_XC,      _YC - _R, 0, lc_cyl)  # south

    # ── Cylinder arcs (CCW: E → N → W → S → E) ───────────────────────────────
    a1 = gmsh.model.geo.addCircleArc(pe, pc, pn)  # east  → north
    a2 = gmsh.model.geo.addCircleArc(pn, pc, pw)  # north → west
    a3 = gmsh.model.geo.addCircleArc(pw, pc, ps)  # west  → south
    a4 = gmsh.model.geo.addCircleArc(ps, pc, pe)  # south → east

    # ── Curve loops and surface ───────────────────────────────────────────────
    cl_outer = gmsh.model.geo.addCurveLoop([l_bottom, l_outlet, l_top, l_inlet])
    cl_cyl   = gmsh.model.geo.addCurveLoop([a1, a2, a3, a4])
    surf     = gmsh.model.geo.addPlaneSurface([cl_outer, cl_cyl])

    gmsh.model.geo.synchronize()
    return surf, (l_inlet, l_outlet, l_bottom, l_top, a1, a2, a3, a4)


def set_refinement_field(arcs: tuple) -> None:
    """
    Add a Distance+Threshold field for smooth size transition near the cylinder.

    Size varies from lc_cyl at the cylinder surface to lc_far at 0.15 m distance
    (approximately 3 cylinder radii), which adequately resolves the wake region.
    """
    a1, a2, a3, a4 = arcs[4], arcs[5], arcs[6], arcs[7]  # cylinder arc tags

    dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist, "CurvesList", [a1, a2, a3, a4])
    gmsh.model.mesh.field.setNumber(dist,  "Sampling",   200)

    thresh = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(thresh, "InField",  dist)
    gmsh.model.mesh.field.setNumber(thresh, "SizeMin",  0.01)   # lc_cyl applied here
    gmsh.model.mesh.field.setNumber(thresh, "SizeMax",  0.05)   # lc_far applied here
    gmsh.model.mesh.field.setNumber(thresh, "DistMin",  0.0)
    gmsh.model.mesh.field.setNumber(thresh, "DistMax",  0.15)

    gmsh.model.mesh.field.setAsBackgroundMesh(thresh)


def add_physical_groups(surf: int, tags: tuple) -> None:
    """
    Register physical groups so libMesh GmshIO assigns boundary IDs.

    Tags:
        1 = inlet, 2 = outlet, 3 = walls, 4 = cylinder, 1 = fluid (2D)
    """
    l_inlet, l_outlet, l_bottom, l_top, a1, a2, a3, a4 = tags

    gmsh.model.addPhysicalGroup(1, [l_inlet],              tag=1, name="inlet")
    gmsh.model.addPhysicalGroup(1, [l_outlet],             tag=2, name="outlet")
    gmsh.model.addPhysicalGroup(1, [l_bottom, l_top],      tag=3, name="walls")
    gmsh.model.addPhysicalGroup(1, [a1, a2, a3, a4],       tag=4, name="cylinder")
    gmsh.model.addPhysicalGroup(2, [surf],                 tag=1, name="fluid")


def generate_mesh(lc_far: float, lc_cyl: float, output: str) -> None:
    """
    Orchestrate geometry, meshing, and export.

    Args:
        lc_far:  Global characteristic element size [m].
        lc_cyl:  Element size at the cylinder surface [m].
        output:  Destination path for the .msh file.
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("channel")

    surf, tags = build_geometry(lc_far, lc_cyl)
    set_refinement_field(tags)
    add_physical_groups(surf, tags)

    # Frontal-Delaunay algorithm gives good quality triangles near curved boundaries
    gmsh.option.setNumber("Mesh.Algorithm",    6)
    gmsh.option.setNumber("Mesh.RecombineAll", 0)   # triangles only
    gmsh.model.mesh.generate(2)

    # Export as Gmsh format 2.2 ASCII — required for libMesh GmshIO
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(output)

    # Summary statistics
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(dim=2, tag=-1)
    n_tri = sum(len(t) for t in elem_tags)
    print(f"\nMesh written: {output}")
    print(f"  Nodes:     {len(node_tags)}")
    print(f"  Triangles: {n_tri}")

    gmsh.finalize()


if __name__ == "__main__":
    args = parse_args()
    generate_mesh(lc_far=args.lc_far, lc_cyl=args.lc_cyl, output=args.output)
