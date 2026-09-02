// Project 01 - Hypersonic Blunt Body
// STAGE 4a (revised): Build a thin FLAT slab (straight extrude, not
// rotational) from the simple domain-outline geometry. This is purely a
// container to get the 2D mesh topology into OpenFOAM via gmshToFoam;
// the actual wedge revolution will be done afterward by OpenFOAM's own
// extrudeMesh utility (extrudeModel wedge), not by Gmsh's Recombine,
// since Gmsh's 3D Recombine was found to produce stray tetrahedra/pyramids
// instead of clean prisms (Stage 4a v1 finding).

R_n = 0.05;
theta_c = 15 * Pi/180;
R_b = 0.15;

x_t = R_n * (1 - Sin(theta_c));
r_t = R_n * Cos(theta_c);
L = x_t + (R_b - r_t) / Tan(theta_c);

x_upstream = -6 * R_n;
r_outer    = 15 * R_n;
x_downstream = 3 * L;

lc = 0.02;

n_nose = 20;
For i In {0:n_nose}
  xi = i * x_t / n_nose;
  ri = Sqrt(R_n^2 - (xi - R_n)^2);
  p = newp;
  Point(p) = {xi, ri, 0, lc};
  nose_pts[i] = p;
EndFor

p_nose_tip   = nose_pts[0];
p_tangent    = nose_pts[n_nose];
p_base_outer = newp; Point(p_base_outer) = {L, R_b, 0, lc};
p_base_axis  = newp; Point(p_base_axis)  = {L, 0, 0, lc};
p_wake_axis  = newp; Point(p_wake_axis)  = {x_downstream, 0, 0, lc};
p_outlet_top = newp; Point(p_outlet_top) = {x_downstream, r_outer, 0, lc};
p_farfield_up = newp; Point(p_farfield_up) = {x_upstream, r_outer, 0, lc};
p_axis_up    = newp; Point(p_axis_up)    = {x_upstream, 0, 0, lc};

Line(1) = {p_axis_up, p_nose_tip};       // axis, upstream
Spline(2) = nose_pts[];                   // wall: nose
Line(3) = {p_tangent, p_base_outer};      // wall: cone
Line(4) = {p_base_outer, p_base_axis};    // wall: base
Line(5) = {p_base_axis, p_wake_axis};     // axis, downstream
Line(6) = {p_wake_axis, p_outlet_top};    // outlet
Line(7) = {p_outlet_top, p_farfield_up};  // farfield, outer
Line(8) = {p_farfield_up, p_axis_up};     // farfield, upstream

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};
Recombine Surface{1};

// --- Thin FLAT slab: straight extrude, tiny thickness, single layer ---
slab_thickness = 0.001; // 1 mm, arbitrary and small - purely a container

out[] = Extrude {0, 0, slab_thickness} {
  Surface{1}; Layers{1}; Recombine;
};

// --- Physical groups: required for gmshToFoam to assign meaningful patch names ---
Physical Surface("frontPatch") = {1};        // the ORIGINAL surface (z=0)
Physical Surface("backPatch") = {out[0]};    // the extruded-to surface (z=slab_thickness)

// out[] indices for a Surface extrusion: out[0]=top surface, out[1]=volume,
// out[2..]=side surfaces in curve-loop order (1,2,3,4,5,6,7,8)
Physical Surface("wall_nose")     = {out[3]};   // side from curve 2
Physical Surface("wall_cone")     = {out[4]};   // side from curve 3
Physical Surface("wall_base")     = {out[5]};   // side from curve 4
Physical Surface("axis_upstream") = {out[2]};   // side from curve 1
Physical Surface("axis_downstream")={out[6]};   // side from curve 5
Physical Surface("outlet")        = {out[7]};   // side from curve 6
Physical Surface("farfield_outer")= {out[8]};   // side from curve 7
Physical Surface("farfield_upstream")={out[9]}; // side from curve 8

Physical Volume("internal") = {out[1]};

// ... (all content identical up to and including the Physical Volume line) ...

// STAGE 4a FIX: disable Gmsh mesh optimization, which was found to
// perturb frontPatch point z-coordinates by floating-point noise during
// 3D tet-mesh smoothing, breaking exact planarity required by OpenFOAM's
// wedge extrudeModel (confirmed via "Wedge patch not planar" errors and
// a subsequent floating-point-exception crash in extrudeMesh).
Mesh.Optimize = 0;
Mesh.OptimizeNetgen = 0;
Mesh.Smoothing = 0;
