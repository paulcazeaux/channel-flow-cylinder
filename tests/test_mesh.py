#!/usr/bin/env python3
"""
test_mesh.py
Unit tests for Phase 1: Gmsh mesh generation.

Runs generate_mesh.py on a coarse mesh (lc_far=0.1, lc_cyl=0.03) so the
test completes quickly on a workstation, then checks structural validity of
the output.  Run with:
    python -m pytest tests/test_mesh.py -v
or:
    python tests/test_mesh.py
"""

import os
import sys
import tempfile
import unittest

# Locate generate_mesh.py relative to this file
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_MESH_DIR = os.path.join(_ROOT, "meshes")
sys.path.insert(0, _MESH_DIR)

try:
    import gmsh
    _HAS_GMSH = True
except ImportError:
    _HAS_GMSH = False

try:
    from generate_mesh import generate_mesh
    _HAS_GENERATOR = True
except ImportError:
    _HAS_GENERATOR = False


# ── Coarse parameters for fast workstation tests ──────────────────────────────
_LC_FAR  = 0.10
_LC_CYL  = 0.03

# Expected ranges for the coarse mesh (empirically derived)
_MIN_NODES = 500
_MAX_NODES = 8_000
_MIN_TRIS  = 900
_MAX_TRIS  = 16_000

# Minimum acceptable Gmsh element quality (0 = degenerate, 1 = equilateral)
_MIN_QUALITY = 0.05


@unittest.skipUnless(_HAS_GMSH,      "gmsh Python package not available")
@unittest.skipUnless(_HAS_GENERATOR, "generate_mesh module not importable")
class TestMeshGeneration(unittest.TestCase):
    """Validate the channel-cylinder mesh produced by generate_mesh.py."""

    @classmethod
    def setUpClass(cls):
        """Generate the mesh once; all test methods share it."""
        cls.tmp = tempfile.NamedTemporaryFile(suffix=".msh", delete=False)
        cls.msh_path = cls.tmp.name
        cls.tmp.close()
        generate_mesh(lc_far=_LC_FAR, lc_cyl=_LC_CYL, output=cls.msh_path)

        # Re-open mesh in Gmsh for inspection
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(cls.msh_path)

    @classmethod
    def tearDownClass(cls):
        gmsh.finalize()
        os.unlink(cls.msh_path)

    # ── File existence ────────────────────────────────────────────────────────

    def test_file_exists(self):
        """Output .msh file must be created on disk."""
        self.assertTrue(os.path.isfile(self.msh_path),
                        f"Mesh file not found: {self.msh_path}")

    def test_file_nonempty(self):
        """Output file must not be empty."""
        self.assertGreater(os.path.getsize(self.msh_path), 0)

    # ── Node and element counts ───────────────────────────────────────────────

    def test_node_count_in_range(self):
        """Node count must be within expected bounds for the coarse mesh."""
        node_tags, _, _ = gmsh.model.mesh.getNodes()
        n = len(node_tags)
        self.assertGreaterEqual(n, _MIN_NODES, f"Too few nodes: {n}")
        self.assertLessEqual(n,   _MAX_NODES,  f"Too many nodes: {n}")

    def test_triangle_count_in_range(self):
        """Triangle count must be within expected bounds for the coarse mesh."""
        _, elem_tags, _ = gmsh.model.mesh.getElements(dim=2, tag=-1)
        n = sum(len(t) for t in elem_tags)
        self.assertGreaterEqual(n, _MIN_TRIS, f"Too few triangles: {n}")
        self.assertLessEqual(n,   _MAX_TRIS,  f"Too many triangles: {n}")

    # ── Physical boundary tags ────────────────────────────────────────────────

    def _boundary_element_count(self, tag: int) -> int:
        """Return number of 1D boundary elements on a physical group tag."""
        groups = gmsh.model.getPhysicalGroups(dim=1)
        for dim, t in groups:
            if t == tag:
                entities = gmsh.model.getEntitiesForPhysicalGroup(dim, t)
                total = 0
                for e in entities:
                    _, etags, _ = gmsh.model.mesh.getElements(dim=1, tag=e)
                    total += sum(len(x) for x in etags)
                return total
        return 0

    def test_inlet_tag_present(self):
        """Physical tag 1 (inlet) must contain boundary elements."""
        self.assertGreater(self._boundary_element_count(1), 0,
                           "No elements found on inlet (tag 1)")

    def test_outlet_tag_present(self):
        """Physical tag 2 (outlet) must contain boundary elements."""
        self.assertGreater(self._boundary_element_count(2), 0,
                           "No elements found on outlet (tag 2)")

    def test_walls_tag_present(self):
        """Physical tag 3 (walls) must contain boundary elements."""
        self.assertGreater(self._boundary_element_count(3), 0,
                           "No elements found on walls (tag 3)")

    def test_cylinder_tag_present(self):
        """Physical tag 4 (cylinder) must contain boundary elements."""
        self.assertGreater(self._boundary_element_count(4), 0,
                           "No elements found on cylinder (tag 4)")

    # ── Element quality ───────────────────────────────────────────────────────

    def test_no_degenerate_elements(self):
        """Minimum element quality (Gmsh SICN metric) must exceed threshold."""
        entities = gmsh.model.getEntities(dim=2)
        min_q = 1.0
        for _, tag in entities:
            etypes, etags, _ = gmsh.model.mesh.getElements(dim=2, tag=tag)
            if not etypes:
                continue
            for tags in etags:
                q = gmsh.model.mesh.getElementQualities(tags, "minSICN")
                if len(q):
                    min_q = min(min_q, min(q))
        self.assertGreater(min_q, _MIN_QUALITY,
                           f"Degenerate elements found (min quality = {min_q:.4f})")


if __name__ == "__main__":
    unittest.main(verbosity=2)
