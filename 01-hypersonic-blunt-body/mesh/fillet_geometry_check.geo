// Project 01 - Hypersonic Blunt Body
// Cone/base corner fillet - GEOMETRY-ONLY verification (no mesh fields).
// Verifies: dimensions, tangency/continuity, unchanged nose/cone profile,
// watertight closed loop, correct patch identification - BEFORE any
// mesh generation, per explicit instruction.

R_n = 0.05;
theta_c = 15 * Pi/180;
R_b = 0.15;
r_fillet = 0.003;  // explicit geometry parameter, not hard-coded

x_t = R_n * (1 - Sin(theta_c));
r_t = R_n * Cos(theta_c);
L = x_t + (R_b - r_t) / Tan(theta_c);

x_upstream = -6 * R_n;
r_outer    = 15 * R_n;
x_downstream = 3 * L;

lc = 0.02;
lc_fillet_check = 0.0002; // finer, for honest visual verification only

// --- Nose spline: UNCHANGED from baseline geometry ---
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

// --- Fillet geometry derivation ---
phi = Pi/2 - theta_c;
t_len = r_fillet / Tan(phi/2);
d_len = r_fillet / Sin(phi/2);
bisector_angle = 5*Pi/4 + theta_c/2;

T1_x = L - t_len * Cos(theta_c);
T1_r = R_b - t_len * Sin(theta_c);
T2_x = L;
T2_r = R_b - t_len;
O_x  = L + d_len * Cos(bisector_angle);
O_r  = R_b + d_len * Sin(bisector_angle);

Printf("Fillet derivation check:");
Printf("  phi (interior angle) = %g deg", phi * 180/Pi);
Printf("  t_len = %g, d_len = %g", t_len, d_len);
Printf("  T1 = (%g, %g)  [on cone line]", T1_x, T1_r);
Printf("  T2 = (%g, %g)  [on base line]", T2_x, T2_r);
Printf("  O  = (%g, %g)  [fillet centre]", O_x, O_r);
Printf("  dist(O,T1) = %g  (should equal r_fillet = %g)", Sqrt((O_x-T1_x)^2+(O_r-T1_r)^2), r_fillet);
Printf("  dist(O,T2) = %g  (should equal r_fillet = %g)", Sqrt((O_x-T2_x)^2+(O_r-T2_r)^2), r_fillet);

// --- Points ---
p_T1 = newp; Point(p_T1) = {T1_x, T1_r, 0, lc_fillet_check};
p_T2 = newp; Point(p_T2) = {T2_x, T2_r, 0, lc_fillet_check};
p_O  = newp; Point(p_O)  = {O_x, O_r, 0, lc};   // fillet centre (construction point)
p_base_axis  = newp; Point(p_base_axis)  = {L, 0, 0, lc};
p_wake_axis  = newp; Point(p_wake_axis)  = {x_downstream, 0, 0, lc};
p_outlet_top = newp; Point(p_outlet_top) = {x_downstream, r_outer, 0, lc};
p_farfield_up = newp; Point(p_farfield_up) = {x_upstream, r_outer, 0, lc};
p_axis_up    = newp; Point(p_axis_up)    = {x_upstream, 0, 0, lc};

// --- Curves ---
Line(1) = {p_axis_up, p_nose_tip};        // axis, upstream
Spline(2) = nose_pts[];                    // wall: nose (UNCHANGED)
Line(3) = {p_tangent, p_T1};               // wall: cone (shortened, ends at T1)
Circle(9) = {p_T1, p_O, p_T2};             // wall: fillet (NEW)
Line(4) = {p_T2, p_base_axis};             // wall: base (shortened, starts at T2)
Line(5) = {p_base_axis, p_wake_axis};      // axis, downstream
Line(6) = {p_wake_axis, p_outlet_top};     // outlet
Line(7) = {p_outlet_top, p_farfield_up};   // farfield, outer
Line(8) = {p_farfield_up, p_axis_up};      // farfield, upstream

Curve Loop(1) = {1, 2, 3, 9, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

// --- Physical groups (curves, for this geometry-check stage) ---
Physical Line("wall_nose")  = {2};
Physical Line("wall_cone")  = {3};
Physical Line("wall_fillet")= {9};
Physical Line("wall_base")  = {4};
Physical Surface("domain") = {1};
