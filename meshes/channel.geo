// channel.geo
// Parametric geometry for 2D channel flow past a circular cylinder.
// Schafer-Turek (1996) DFG 2D-1 benchmark domain.
//
// Can be meshed directly: gmsh channel.geo -2 -format msh2 -o channel.msh
// Or used via generate_mesh.py for scripted/automated mesh generation.
//
// Physical boundary IDs (match libMesh boundary_id in main.cpp):
//   1 = inlet    (x = 0)
//   2 = outlet   (x = L)
//   3 = walls    (y = 0 and y = H)
//   4 = cylinder surface

// ── Domain parameters ────────────────────────────────────────────────────────
H     = 0.41;    // channel height [m]
L     = 2.2;     // channel length [m]
xc    = 0.2;     // cylinder center x [m]
yc    = 0.2;     // cylinder center y [m]
r     = 0.05;    // cylinder radius [m]

lc_far = 0.05;   // global characteristic element size [m]
lc_cyl = 0.01;   // near-cylinder element size [m]

// ── Channel corners ──────────────────────────────────────────────────────────
Point(1) = {0, 0, 0, lc_far};   // bottom-left
Point(2) = {L, 0, 0, lc_far};   // bottom-right
Point(3) = {L, H, 0, lc_far};   // top-right
Point(4) = {0, H, 0, lc_far};   // top-left

// ── Channel boundary edges ───────────────────────────────────────────────────
Line(1) = {1, 2};   // bottom wall   (→)
Line(2) = {2, 3};   // outlet        (↑)
Line(3) = {3, 4};   // top wall      (←)
Line(4) = {4, 1};   // inlet         (↓)

// ── Cylinder: center (arc pivot) + 4 boundary points ─────────────────────────
Point(5) = {xc,   yc,   0, lc_cyl};   // center (not a mesh boundary node)
Point(6) = {xc+r, yc,   0, lc_cyl};   // east
Point(7) = {xc,   yc+r, 0, lc_cyl};   // north
Point(8) = {xc-r, yc,   0, lc_cyl};   // west
Point(9) = {xc,   yc-r, 0, lc_cyl};   // south

// ── Cylinder arcs (CCW: E → N → W → S → E) ──────────────────────────────────
Circle(5) = {6, 5, 7};   // east  → north
Circle(6) = {7, 5, 8};   // north → west
Circle(7) = {8, 5, 9};   // west  → south
Circle(8) = {9, 5, 6};   // south → east

// ── Curve loops and fluid surface ────────────────────────────────────────────
Curve Loop(1) = {1, 2, 3, 4};      // outer channel boundary (CCW)
Curve Loop(2) = {5, 6, 7, 8};      // cylinder hole (CCW — treated as interior hole)
Plane Surface(1) = {1, 2};          // fluid domain = channel minus cylinder

// ── Mesh refinement: smooth size transition from cylinder to far field ────────
Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8};
Field[1].Sampling   = 200;          // points sampled along each curve

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = lc_cyl;         // element size at cylinder surface
Field[2].SizeMax  = lc_far;         // element size in the far field
Field[2].DistMin  = 0.0;            // start coarsening immediately outside cylinder
Field[2].DistMax  = 0.15;           // far-field size reached at 3 radii distance

Background Field = 2;

// ── Physical groups ───────────────────────────────────────────────────────────
Physical Curve("inlet",    1) = {4};
Physical Curve("outlet",   2) = {2};
Physical Curve("walls",    3) = {1, 3};
Physical Curve("cylinder", 4) = {5, 6, 7, 8};
Physical Surface("fluid",  1) = {1};

// ── Mesh algorithm options ────────────────────────────────────────────────────
Mesh.Algorithm   = 6;   // Frontal-Delaunay: good quality near curved boundaries
Mesh.RecombineAll = 0;  // Keep triangles (required for Taylor-Hood P2/P1 in libMesh)
