// Project 01 - Hypersonic Blunt Body
// STAGE 4a v2: Gmsh native rotational extrude, symmetric wedge built by
// (1) rotating the 2D surface to theta=-1deg first, then (2) extruding
// via Rotate by the full +2deg to sweep it to theta=+1deg. This sequencing
// is a documented working pattern (CFD-Online, Eric Tridas) for producing
// OpenFOAM-valid wedge meshes directly via gmshToFoam, without needing
// OpenFOAM's own extrudeMesh utility (which hit repeated, unresolved
// face-orientation crashes in Stage 4a v1-alt attempts).

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

Line(1) = {p_axis_up, p_nose_tip};
Spline(2) = nose_pts[];
Line(3) = {p_tangent, p_base_outer};
Line(4) = {p_base_outer, p_base_axis};
Line(5) = {p_base_axis, p_wake_axis};
Line(6) = {p_wake_axis, p_outlet_top};
Line(7) = {p_outlet_top, p_farfield_up};
Line(8) = {p_farfield_up, p_axis_up};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};
Recombine Surface{1};

half_angle = 1 * Pi/180;
full_angle = 2 * half_angle;

// Step 1: rotate the ENTIRE surface (all its points) to theta = -half_angle
Rotate { {1,0,0}, {0,0,0}, -half_angle } { Surface{1}; }

// Step 2: rotationally extrude by the FULL angle, sweeping from
// theta=-half_angle to theta=+half_angle
out[] = Extrude { {1,0,0}, {0,0,0}, full_angle } {
  Surface{1}; Layers{1}; Recombine;
};

Coherence;   // merge coincident/duplicate points (e.g. axis points that
             // map to the same location under both +/-1deg rotation,
             // since r=0 points are invariant under rotation about r=0)

// --- Physical groups ---
Physical Surface("frontWedge") = {1};      // now sitting at theta = -1deg
Physical Surface("backWedge")  = {out[0]}; // now sitting at theta = +1deg

// Curves 1 (axis_upstream) and 5 (axis_downstream) lie ON the rotation
// axis and produce NO swept surface under rotational extrude (this is
// correct axisymmetric geometry, not an error) - so we do not tag them
// here. Their faces will appear in gmshToFoam's "defaultFaces" and will
// be retyped to "empty" afterward via createPatch, matching OpenFOAM's
// documented axisymmetric wedge convention (axis -> empty type).
Physical Surface("wall_nose")      = {out[2]};
Physical Surface("wall_cone")      = {out[3]};
Physical Surface("wall_base")      = {out[4]};
Physical Surface("outlet")         = {out[5]};
Physical Surface("farfield_outer") = {out[6]};
Physical Surface("farfield_upstream")={out[7]};

Physical Volume("internal") = {out[1]};

Mesh.Optimize = 0;
Mesh.OptimizeNetgen = 0;
Mesh.Smoothing = 0;
