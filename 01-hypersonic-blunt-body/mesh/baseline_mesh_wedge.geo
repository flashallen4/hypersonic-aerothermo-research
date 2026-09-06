// Project 01 - Hypersonic Blunt Body
// STAGE 3b: Full baseline mesh - boundary layer (wall) + shock-layer
// refinement zone (Box field, Option C) + coarse farfield.
//
// This is now the primary baseline mesh script. All sizing parameters are
// expressed in terms of R_n, theta_c, R_b so this file's structure can be
// reused for the nose-radius sweep (only R_n needs to change).

R_n = 0.05;
theta_c = 15 * Pi/180;
R_b = 0.15;

x_t = R_n * (1 - Sin(theta_c));
r_t = R_n * Cos(theta_c);
L = x_t + (R_b - r_t) / Tan(theta_c);

x_upstream = -6 * R_n;
r_outer    = 15 * R_n;
x_outlet = L + 0.3 * R_n;

lc_far   = 0.03;
lc_wall  = 0.001;
lc_shock = R_n / 12.0;

// --- Body profile points (nose spline) - fine lc_wall ---
n_nose = 40;
For i In {0:n_nose}
  xi = i * x_t / n_nose;
  ri = Sqrt(R_n^2 - (xi - R_n)^2);
  p = newp;
  Point(p) = {xi, ri, 0, lc_wall};
  nose_pts[i] = p;
EndFor

p_nose_tip   = nose_pts[0];
p_tangent    = nose_pts[n_nose];
p_base_outer = newp; Point(p_base_outer) = {L, R_b, 0, lc_wall};
p_base_axis  = newp; Point(p_base_axis)  = {L, 0, 0, lc_far};
p_short_axis = newp; Point(p_short_axis) = {x_outlet, 0, 0, lc_far};
p_outlet_top = newp; Point(p_outlet_top) = {x_outlet, r_outer, 0, lc_far};
p_farfield_up = newp; Point(p_farfield_up) = {x_upstream, r_outer, 0, lc_far};
p_axis_up    = newp; Point(p_axis_up)    = {x_upstream, 0, 0, lc_far};

Line(1) = {p_axis_up, p_nose_tip};
Spline(2) = nose_pts[];
Line(3) = {p_tangent, p_base_outer};
Line(4) = {p_base_outer, p_base_axis};
Line(5) = {p_base_axis, p_short_axis};
Line(6) = {p_short_axis, p_outlet_top};
Line(7) = {p_outlet_top, p_farfield_up};
Line(8) = {p_farfield_up, p_axis_up};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

// --- Field 1: Boundary Layer (wall-normal growth) ---
Field[1] = BoundaryLayer;
Field[1].EdgesList = {2, 3, 4};
Field[1].hwall_n = 5.607e-6;
Field[1].ratio = 1.12;
Field[1].thickness = 0.006;
Field[1].Quads = 1;
Field[1].NbLayers = 40;
Field[1].FanPointsList = {p_nose_tip};
BoundaryLayer Field = 1;

// --- Field 2: Shock-layer refinement box (Option C) ---
Field[2] = Box;
Field[2].XMin = -3 * R_n;
Field[2].XMax = L;
Field[2].YMin = 0;
Field[2].YMax = R_b + 4 * R_n;
Field[2].VIn = lc_shock;
Field[2].VOut = lc_far;
Field[2].Thickness = 2 * R_n;

// --- Background mesh size field = shock-layer box (governs freestream region;
//     BoundaryLayer field separately governs near-wall growth) ---
Background Field = 2;

Mesh.Algorithm = 6;

// ============================================================================
// Controlled experiment: sharp corner + extended EdgesList {2,3,4}.
// Index mapping for 8-curve loop {1,2,3,4,5,6,7,8}: curves 1,5 axis-touching
// (no swept surface). out[2]=curve2(nose), out[3]=curve3(cone),
// out[4]=curve4(base), out[5]=curve6(outlet), out[6]=curve7(farfield_outer),
// out[7]=curve8(farfield_upstream).
// ============================================================================

half_angle = 1 * Pi/180;
full_angle = 2 * half_angle;

Rotate { {1,0,0}, {0,0,0}, -half_angle } { Surface{1}; }

out[] = Extrude { {1,0,0}, {0,0,0}, full_angle } {
  Surface{1}; Layers{1}; Recombine;
};

Physical Surface("frontWedge") = {1};
Physical Surface("backWedge")  = {out[0]};

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
