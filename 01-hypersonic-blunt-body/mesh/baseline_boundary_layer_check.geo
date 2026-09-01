// Project 01 - Hypersonic Blunt Body
// STAGE 3a (revised): Boundary-layer field near the wall, rest of domain coarse.
//
// FIX applied: BoundaryLayer field only controls growth NORMAL to the wall.
// Tangential point spacing along the wall curves is governed by the
// characteristic length assigned to the points themselves during 1D meshing
// (a known Gmsh behavior - BoundaryLayer field has no NumPointsPerCurve
// control, unlike Distance/AttractorAnisoCurve fields). Wall points must
// therefore get a FINE lc directly, separate from the coarse farfield lc.

R_n = 0.05;
theta_c = 15 * Pi/180;
R_b = 0.15;

x_t = R_n * (1 - Sin(theta_c));
r_t = R_n * Cos(theta_c);
L = x_t + (R_b - r_t) / Tan(theta_c);

x_upstream = -6 * R_n;
r_outer    = 15 * R_n;
x_downstream = 3 * L;

lc_far  = 0.03;    // coarse background size, far from wall
lc_wall = 0.001;   // fine TANGENTIAL spacing along wall curves (1mm along-wall resolution;
                    // normal-direction resolution is governed separately by the
                    // BoundaryLayer field's hwall_n/ratio/thickness, not by this value)

// --- Body profile points (nose spline) - FINE lc_wall ---
n_nose = 40; // increased point count for smoother spline with fine lc
For i In {0:n_nose}
  xi = i * x_t / n_nose;
  ri = Sqrt(R_n^2 - (xi - R_n)^2);
  p = newp;
  Point(p) = {xi, ri, 0, lc_wall};
  nose_pts[i] = p;
EndFor

p_nose_tip   = nose_pts[0];
p_tangent    = nose_pts[n_nose];
p_base_outer = newp; Point(p_base_outer) = {L, R_b, 0, lc_wall};      // fine: end of cone wall
p_base_axis  = newp; Point(p_base_axis)  = {L, 0, 0, lc_far};          // farfield-adjacent
p_wake_axis  = newp; Point(p_wake_axis)  = {x_downstream, 0, 0, lc_far};
p_outlet_top = newp; Point(p_outlet_top) = {x_downstream, r_outer, 0, lc_far};
p_farfield_up = newp; Point(p_farfield_up) = {x_upstream, r_outer, 0, lc_far};
p_axis_up    = newp; Point(p_axis_up)    = {x_upstream, 0, 0, lc_far};

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

// --- Boundary Layer field: applied to wall curves only (nose=2, cone=3) ---
Field[1] = BoundaryLayer;
Field[1].EdgesList = {2, 3};
Field[1].hwall_n = 5.607e-6;
Field[1].ratio = 1.12;
Field[1].thickness = 0.006;
Field[1].Quads = 1;
Field[1].NbLayers = 40;
Field[1].FanPointsList = {p_nose_tip};  // fan treatment at the stagnation-point corner
                                          // (curve 1/axis meets curve 2/nose here)
BoundaryLayer Field = 1;

Mesh.Algorithm = 6;
