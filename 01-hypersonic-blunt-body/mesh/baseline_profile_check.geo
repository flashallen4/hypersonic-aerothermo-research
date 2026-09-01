// Project 01 - Hypersonic Blunt Body
// STAGE 1: Body profile geometry check only.
// Purpose: verify the spherically-blunted-cone profile renders correctly
// in Gmsh before adding BL mesh, farfield domain, or wedge extrusion.
//
// Parameters match scripts/geometry/blunted_cone_profile.py baseline case.

// --- Baseline parameters ---
R_n = 0.05;           // nose radius [m]
theta_c = 15 * Pi/180; // cone half-angle [rad]
R_b = 0.15;           // base radius [m]

// --- Derived tangent point and body length (must match Python script) ---
x_t = R_n * (1 - Sin(theta_c));
r_t = R_n * Cos(theta_c);
L = x_t + (R_b - r_t) / Tan(theta_c);

Printf("x_t = %g", x_t);
Printf("r_t = %g", r_t);
Printf("L   = %g", L);

// --- Mesh element size for this geometry-only check (coarse, just visual) ---
lc = 0.005;

// --- Nose cap: discretized points along the sphere, then splined ---
n_nose_pts = 15;
For i In {0:n_nose_pts}
  t = i * (Pi/2) / n_nose_pts;  // parametrize sphere cap from tip (t=0) to tangent (t~ related to theta_c)
  // Sphere center at (R_n, 0); tip at t=90deg from center-to-axis direction... 
  // Use explicit profile formula instead for exact match to Python script:
EndFor

// Simpler and exact: sample x from 0 to x_t, compute r(x) = sqrt(R_n^2 - (x-R_n)^2)
nose_x[] = {};
nose_r[] = {};
n_nose = 20;
For i In {0:n_nose}
  xi = i * x_t / n_nose;
  ri = Sqrt(R_n^2 - (xi - R_n)^2);
  p = newp;
  Point(p) = {xi, ri, 0, lc};
  nose_x[i] = xi;
  nose_r[i] = ri;
  If (i == 0)
    nose_pt_start = p;
  EndIf
  If (i == n_nose)
    nose_pt_end = p;
  EndIf
  nose_pts[i] = p;
EndFor

// Cone line: from tangent point to base
p_cone_end = newp;
Point(p_cone_end) = {L, R_b, 0, lc};

// --- Curves ---
Spline(1) = nose_pts[];
Line(2) = {nose_pt_end, p_cone_end};

// --- Visual check only: no surfaces/mesh yet ---
