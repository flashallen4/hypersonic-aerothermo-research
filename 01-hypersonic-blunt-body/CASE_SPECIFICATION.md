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

| Quantity | Symbol | Baseline Value | Units |
|---|---|---|---|
| Nose radius | R_n | 0.05 | m |
| Cone half-angle | θ_c | 15 | deg |
| Base radius | R_b | 0.15 | m |
| Overall body length | L | computed from R_n, θ_c, R_b | m |

Body length L is a dependent quantity: the spherical nose cap blends tangentially
into the conical afterbody, and the cone extends until it reaches R_b. Exact L
will be computed and reported by the geometry-generation script.

## Nose Radius Sweep (tentative — final set confirmed after baseline is verified/validated)

Baseline case (R_n = 0.05 m) will be built, meshed, verified (grid convergence),
and validated (vs. Fay–Riddell and Billig correlations) FIRST. Only after the
baseline is trustworthy will the sweep be automated. Tentative sweep values:

- R_n = 0.02 m
- R_n = 0.05 m (baseline)
- R_n = 0.10 m
- R_n = 0.20 m

Cone half-angle (15°) and base radius (0.15 m) held fixed across the sweep.

## Wall Boundary Condition

- Fixed isothermal wall temperature: T_w = 300 K

## Domain and Boundary Conditions (to be finalized at meshing stage)

- 2D-axisymmetric wedge domain (OpenFOAM axisymmetric wedge convention)
- Farfield/inlet: freestream supersonic inflow (fixed M∞, T∞, p∞)
- Outlet: supersonic outflow (zero-gradient, since flow is supersonic)
- Body surface: no-slip, fixed-temperature wall (300 K)
- Axis: axisymmetric wedge boundary condition
- Domain extents and farfield distance to be set to avoid boundary influence on
  bow shock — confirmed via a domain-independence check.

## Solver

- `rhoCentralFoam` (OpenFOAM 11), steady-state operation via appropriate
  pseudo-transient/local-time-stepping approach or steady-state analog
  (exact time-stepping/convergence strategy to be finalized at solver-setup
  stage and documented there).

## Status

Baseline parameters confirmed. Geometry generation not yet started.
