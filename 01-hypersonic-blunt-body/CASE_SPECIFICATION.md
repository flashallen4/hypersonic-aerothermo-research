# Project 01 — Baseline Case Specification

## Freestream Conditions (US Standard Atmosphere, 30 km altitude)

| Quantity | Symbol | Value | Units |
|---|---|---|---|
| Altitude | — | 30 | km |
| Static temperature | T∞ | 226.5 | K |
| Static pressure | p∞ | 1197 | Pa |
| Density | ρ∞ | 1.841×10⁻² | kg/m³ |
| Ratio of specific heats | γ | 1.4 | — |
| Specific gas constant (air) | R | 287.0 | J/(kg·K) |
| Speed of sound | a∞ = √(γRT∞) | ≈ 301.7 | m/s |
| Freestream Mach number | M∞ | 7.0 | — |
| Freestream velocity | U∞ = M∞·a∞ | ≈ 2112 | m/s |

Dynamic viscosity (Sutherland's law, for laminar Re/boundary-layer resolution
during meshing) to be computed in the meshing stage from T∞.

## Baseline Geometry — Spherically-Blunted Cone (Axisymmetric)

Profile: spherical nose cap (radius R_n) blended tangentially into a conical
afterbody (half-angle θ_c), extending to base radius R_b.

**Governing relations:**

Spherical cap: r(x) = sqrt(R_n^2 - (x - R_n)^2), 0 <= x <= x_t
Tangent point: r_t = R_n * cos(theta_c)
x_t = R_n * (1 - sin(theta_c))
Conical afterbody: r(x) = r_t + (x - x_t) * tan(theta_c), x_t <= x <= L
Body length: L = x_t + (R_b - r_t) / tan(theta_c)


| Quantity | Symbol | Baseline Value | Units |
|---|---|---|---|
| Nose radius | R_n | 0.05 | m |
| Cone half-angle | θ_c | 15 | deg |
| Base radius | R_b | 0.15 | m |
| Tangent point (axial) | x_t | 0.03706 | m |
| Tangent point (radial) | r_t | 0.04830 | m |
| Overall body length | L | 0.4166 | m |

## Nose Radius Sweep

Baseline case (R_n = 0.05 m) will be built, meshed, verified (grid convergence),
and validated (vs. Fay–Riddell and Billig correlations) FIRST. Only after the
baseline is trustworthy will the sweep be automated.

Cone half-angle (15°) and base radius (0.15 m) held fixed across all sweep cases
(single-variable study).

| R_n (m) | x_t (m) | r_t (m) | L (m) | Notes |
|---|---|---|---|---|
| 0.02 | 0.01482 | 0.01932 | 0.50253 | |
| 0.05 | 0.03706 | 0.04830 | 0.41662 | Baseline |
| 0.10 | 0.07412 | 0.09659 | 0.27344 | |
| 0.14 | 0.10377 | 0.13523 | 0.15889 | Bluntest case; short (~5.5 cm) conical afterbody |

Values in this table are generated and verified by
`scripts/geometry/blunted_cone_profile.py` (single source of truth for all
profile geometry). Earlier hand-calculated values for R_n=0.02 and R_n=0.10
contained arithmetic errors, caught by cross-checking against this script.

**Note on sweep upper bound:** R_n = 0.20 m was initially proposed but found
geometrically invalid at θ_c = 15°, R_b = 0.15 m — the tangent-point radius
r_t = 0.1932 m already exceeds R_b before the sphere reaches tangency with the
cone, meaning no valid conical afterbody exists at that combination. Sweep
capped at R_n = 0.14 m to preserve a strictly single-variable (nose-radius-only)
study with θ_c and R_b held fixed for all cases.

## Wall Boundary Condition

- Fixed isothermal wall temperature: T_w = 300 K

## Domain, Mesh Sizing, and Boundary Conditions

- 2D-axisymmetric wedge domain (OpenFOAM axisymmetric wedge convention,
  wedge half-angle ~2-5 deg)
- Farfield/inlet: freestream supersonic inflow (fixed M∞, T∞, p∞)
- Outlet: supersonic outflow (zero-gradient, since flow is supersonic)
- Body surface: no-slip, fixed-temperature wall (300 K)
- Axis: axisymmetric wedge boundary condition

**Domain extents (baseline case, R_n = 0.05 m):**

| Boundary | Distance | Rationale |
|---|---|---|
| Upstream farfield | 6 × R_n = 0.30 m ahead of nose | Comfortable margin ahead of bow shock stand-off |
| Outer radial farfield | 15 × R_n = 0.75 m from axis | Avoid confinement of shock layer / outer flow |
| Downstream extent | 2 × L past base (~3 × L total) | Allow wake/expansion to develop away from outlet |

To be confirmed via a domain-independence check (part of verification stage):
if results change meaningfully when farfield boundaries are pushed further out,
the domain will be enlarged.

**Near-wall mesh sizing (baseline case):**

Computed via `scripts/mesh/boundary_layer_estimate.py`, using Sutherland's law
for viscosity and the Eckert reference-temperature method for compressible
laminar boundary-layer thickness estimation (a standard first-order sizing
approximation — NOT used as a validation or physics result; the CFD solution
itself is what will be trusted for actual boundary-layer physics).

| Quantity | Value |
|---|---|
| Adiabatic wall temperature (recovery) | 2113.24 K |
| Eckert reference temperature | 678.33 K |
| Estimated BL thickness near stagnation (thinnest region) | 0.561 mm |
| Estimated BL thickness near base (thickest region) | 5.015 mm |
| Recommended first-cell wall-normal height | 5.61 μm (5.607×10⁻⁶ m) |
| Wall-normal geometric growth ratio | 1.12 |
| Cells spanning thinnest BL region | ~23 |

This first-cell height and growth ratio apply to the baseline case (R_n = 0.05 m).
Sweep cases with different R_n will need this recomputed, since stagnation-region
BL thickness scales with local geometry — this will be automated once the
baseline mesh/solver setup is verified.

## Solver

- `rhoCentralFoam` (OpenFOAM 11), steady-state operation via appropriate
  pseudo-transient/local-time-stepping approach or steady-state analog
  (exact time-stepping/convergence strategy to be finalized at solver-setup
  stage and documented there).

## Status

Baseline and sweep geometry parameters confirmed (sweep upper bound corrected
from 0.20 m to 0.14 m after tangency-validity check). Geometry generation not
yet started.
