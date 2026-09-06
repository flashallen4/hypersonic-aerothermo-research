// Minimal test: two BoundaryLayer fields with DIFFERENT growth parameters,
// combined via a Min field, applied to two different edges of a simple
// rectangle. Purpose: verify Gmsh accepts this combination and produces
// visibly different near-wall spacing on each edge, before applying the
// same technique to the real Project 01 baseline mesh.

lc_far = 0.05;
lc_wall = 0.002;   // fine tangential spacing on BL walls, mirrors real mesh's lc_wall pattern

// Rectangle: (0,0) to (1, 0.5)
// p1, p2 define bottom wall (fine BL) -> lc_wall
// p1, p4 define left wall (coarse BL) -> lc_wall
// p3 is the far corner, no BL -> lc_far
p1 = newp; Point(p1) = {0, 0, 0, lc_wall};
p2 = newp; Point(p2) = {1, 0, 0, lc_wall};
p3 = newp; Point(p3) = {1, 0.5, 0, lc_far};
p4 = newp; Point(p4) = {0, 0.5, 0, lc_wall};

l1 = newl; Line(l1) = {p1, p2};  // bottom wall - FINE BL (mimics wall_cone)
l2 = newl; Line(l2) = {p2, p3};  // right - no BL
l3 = newl; Line(l3) = {p3, p4};  // top - no BL
l4 = newl; Line(l4) = {p4, p1};  // left wall - COARSE BL (mimics wall_base)

Curve Loop(1) = {l1, l2, l3, l4};
Plane Surface(1) = {1};

// --- BoundaryLayer field 1: FINE, on bottom wall (l1) ---
Field[1] = BoundaryLayer;
Field[1].EdgesList = {l1};
Field[1].hwall_n = 5.607e-6;
Field[1].ratio = 1.12;
Field[1].thickness = 0.006;
Field[1].hfar = lc_far;
Field[1].Quads = 1;
Field[1].NbLayers = 40;

// --- BoundaryLayer field 2: COARSE, on left wall (l4) ---
Field[2] = BoundaryLayer;
Field[2].EdgesList = {l4};
Field[2].hwall_n = 5.607e-5;   // 10x coarser first-cell height
Field[2].ratio = 1.2;          // faster growth
Field[2].thickness = 0.006;
Field[2].hfar = lc_far;
Field[2].Quads = 1;
Field[2].NbLayers = 20;

// --- Combine via Min field, per documented working pattern ---
Field[3] = Min;
Field[3].FieldsList = {1, 2};

Background Field = 3;
