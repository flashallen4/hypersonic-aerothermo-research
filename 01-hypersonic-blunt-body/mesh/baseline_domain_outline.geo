// Project 01 - Hypersonic Blunt Body
// STAGE 2: Full domain outline (body + farfield boundary), no BL/sizing yet.
// Verifies the closed-loop domain topology before adding mesh density control.

// --- Baseline body parameters (must match blunted_cone_profile.py) ---
R_n = 0.05;
theta_c = 15 * Pi/180;
R_b = 0.15;

x_t = R_n * (1 - Sin(theta_c));
r_t = R_n * Cos(theta_c);
L = x_t + (R_b - r_t) / Tan(theta_c);

// --- Domain sizing (per CASE_SPECIFICATION.md) ---
x_upstream = -6 * R_n;      // -0.30 m
r_outer    = 15 * R_n;      // 0.75 m
x_downstream = 3 * L;       // L + 2L past base

Printf("x_t=%g r_t=%g L=%g", x_t, r_t, L);
Printf("x_upstream=%g r_outer=%g x_downstream=%g", x_upstream, r_outer, x_downstream);

lc = 0.02; // coarse, outline-check only

// --- Body profile points (nose spline) ---
n_nose = 20;
For i In {0:n_nose}
  xi = i * x_t / n_nose;
  ri = Sqrt(R_n^2 - (xi - R_n)^2);
  p = newp;
  Point(p) = {xi, ri, 0, lc};
  nose_pts[i] = p;
EndFor

// --- Key corner points ---
p_nose_tip   = nose_pts[0];             // (0, 0)
p_tangent    = nose_pts[n_nose];        // tangent point
p_base_outer = newp; Point(p_base_outer) = {L, R_b, 0, lc};        // (L, R_b)
p_base_axis  = newp; Point(p_base_axis)  = {L, 0, 0, lc};          // (L, 0)
p_wake_axis  = newp; Point(p_wake_axis)  = {x_downstream, 0, 0, lc}; // (3L, 0)
p_outlet_top = newp; Point(p_outlet_top) = {x_downstream, r_outer, 0, lc}; // (3L, r_outer)
p_farfield_up = newp; Point(p_farfield_up) = {x_upstream, r_outer, 0, lc}; // (-0.30, r_outer)
p_axis_up    = newp; Point(p_axis_up)    = {x_upstream, 0, 0, lc}; // (-0.30, 0)

// --- Curves (closed loop, counterclockwise) ---
Line(1) = {p_axis_up, p_nose_tip};          // axis, upstream
Spline(2) = nose_pts[];                      // wall: nose cap
Line(3) = {p_tangent, p_base_outer};         // wall: cone
Line(4) = {p_base_outer, p_base_axis};       // wall: flat base
Line(5) = {p_base_axis, p_wake_axis};        // axis, downstream/wake
Line(6) = {p_wake_axis, p_outlet_top};       // outlet
Line(7) = {p_outlet_top, p_farfield_up};     // farfield, outer radial
Line(8) = {p_farfield_up, p_axis_up};        // farfield, upstream

// --- Closed loop + surface (to confirm topology is valid/closed) ---
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};
