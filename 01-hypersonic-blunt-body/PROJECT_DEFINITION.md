# Project 01 — Hypersonic Blunt-Body Aerodynamics and Aerothermal Characterization

## Research Area
Hypersonic aerodynamics / aerothermodynamics

## Engineering Problem

Hypersonic vehicles with blunt forebodies (re-entry capsules, blunt leading edges,
interceptor noses) generate a strong detached bow shock. This shock layer governs
surface pressure distribution, drag, and aerodynamic heating. Nose bluntness is a
primary design lever for managing heat flux, but the relationship between nose
radius, shock stand-off distance, drag, and heating is nonlinear and must be
quantified computationally for a given geometry family and flow regime.

## Research Question

How does nose (bluntness) radius influence shock stand-off distance, surface
pressure distribution, drag, and the distribution/peak of surface heat flux on an
axisymmetric hypersonic blunt body at fixed freestream conditions?

## Engineering Objective

Establish a verified, mesh-independent, and correlation-validated compressible CFD
methodology in OpenFOAM 11 for hypersonic blunt-body flow, then apply it to quantify
the nose-radius → aerothermal-loading trade-off via a single-variable parametric
sweep.

## Hypothesis

Increasing nose radius (at fixed cone half-angle and base radius) will:
- increase shock stand-off distance,
- decrease peak stagnation-point heat flux, consistent with Fay–Riddell scaling
  (q̇_stag ∝ R⁻¹ᐟ²),
- increase pressure drag due to increased frontal bluntness.

These trends are hypotheses to be tested against our own CFD results, not assumed
conclusions.

## Physical Regime and Major Modelling Decisions

| Decision | Choice | Justification |
|---|---|---|
| Freestream regime | Mach ≈ 7, perfect gas (γ = 1.4), ~30 km altitude | Genuinely hypersonic (M > 5) while remaining thermally/chemically frozen — tractable without real-gas chemistry, which OpenFOAM 11 does not natively support |
| Geometry | Spherically-blunted cone, axisymmetric | Reusable methodology for Project 10 (aerothermal optimization); still has established correlations for V&V |
| Solver | `rhoCentralFoam` | Density-based, shock-capturing — standard OpenFOAM solver for compressible supersonic/hypersonic flow |
| Turbulence | Laminar | Keeps scope honest; laminar stagnation heating is exactly the regime Fay–Riddell was derived for; turbulent/transitional heating is Project 03's scope |
| Dimensionality | 2D-axisymmetric (wedge mesh) | Physically appropriate at zero yaw; research question does not involve angle-of-attack effects |
| Wall thermal BC | Fixed isothermal wall, 300 K | Isolates the aerodynamic effect of bluntness on heating; avoids conflating with wall-temperature feedback (reserved for Project 04, TPS) |
| Time treatment | Steady-state | Blunt-body hypersonic flow reaches a steady shock structure under fixed freestream conditions |
| Radiation | Not modeled | Documented limitation |
| Ablation | Not modeled — rigid, non-ablating wall | Documented limitation |

## Major Design Variable (Parametric Sweep)

- **Primary:** Nose (bluntness) radius, R
- **Held fixed during Project 01 sweep:** cone half-angle (15°), base radius (0.15 m)

## Constraints

- Steady, laminar, calorically perfect gas, no radiation, no ablation, no
  real-gas chemistry, fixed isothermal wall.
- 2D-axisymmetric domain only (no 3D / angle-of-attack cases in this project).

## Required Outputs

1. Shock stand-off distance δ as a function of nose radius R
2. Surface pressure coefficient distribution Cp(s)
3. Surface heat flux distribution q̇(s), with emphasis on stagnation value q̇_stag
4. Drag coefficient Cd(R)
5. Validation comparison: q̇_stag(R) vs. Fay–Riddell correlation
6. Verification/validation comparison: δ(R) vs. Billig's empirical correlation

## Verification Strategy

- Grid convergence study on the baseline case: ≥3 mesh densities, tracking
  stagnation-point heat flux and shock stand-off distance. Richardson
  extrapolation / GCI reported where practical.
- Global mass conservation check (inflow vs. outflow mass flux).

## Validation Strategy

- Stagnation-point heat flux compared against the **Fay–Riddell** laminar
  stagnation-point heating correlation.
- Shock stand-off distance compared against **Billig's** empirical correlation
  for blunt-body shock stand-off.
- This is validation against established engineering correlations, not against
  a specific experimental dataset. If a suitable digitized experimental dataset
  (e.g., Lobb 1964) is later located and verified as genuine, it may be added;
  none is assumed or fabricated here.

## Limitations (documented up front, to be revisited at project completion)

- Perfect-gas (calorically perfect, γ = 1.4) assumption — no real-gas
  dissociation/ionization effects relevant at true re-entry Mach numbers.
- No radiative heat transfer (surface-to-surroundings or shock-layer radiation).
- No ablation or surface mass loss.
- Laminar flow assumption — no boundary-layer transition or turbulent heating
  augmentation.
- Fixed isothermal wall — does not capture radiative-equilibrium wall-temperature
  feedback relevant to real TPS design (deferred to Project 04).
- 2D-axisymmetric only — no angle-of-attack or 3D asymmetric effects.
- Validation is against analytical/empirical correlations, not direct
  experimental data.

## Status

**In progress** — specification approved, directory structure created.
Geometry, mesh, and solver setup not yet started.
