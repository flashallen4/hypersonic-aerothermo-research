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

---

## UPDATE — 2026-09-06: Corner Topology Fix, Domain Restructuring, and Corner Instability Investigation

### Summary
Following the geometry/mesh-sizing work documented above, this session investigated and resolved a mesh-topology defect at the cone/base corner, adopted a revised baseline mesh, tested a domain-shortening change to address a separate numerical instability, and conducted a controlled investigation into that instability's root cause. The instability itself remains **unresolved and is accepted as a known, documented limitation** — see CASE_SPECIFICATION.md for full technical detail.

### Key decisions made this session
1. **Baseline mesh adopted: "Case 4"** — original sharp-corner geometry (no fillet) with the BoundaryLayer field's `EdgesList` extended from `{2,3}` to `{2,3,4}` (including `wall_base`). This was established via a rigorous 2×2 controlled experiment as the cause of a previously-observed open-cell/misoriented-face defect at the cone/base corner, independent of wall geometry. A 3mm fillet, tested in parallel, was confirmed **not required** to resolve this defect.
2. **Downstream domain extent shortened**: outlet moved from 3×L (~1.25m total domain) to L+0.3×R_n (~1.5cm past the base), to address a separate negative-pressure/Courant-instability found during solver-viability testing on the clean Case 4 mesh. This required two failed geometric attempts (collinear-edge and near-parallel-edge degeneracies) before a working right-angle-corner solution was found.
3. **A distinct instability at the cone/base corner region was identified, investigated, and NOT resolved.** It is documented as a known limitation (see below and CASE_SPECIFICATION.md). Domain-shortening was conclusively shown NOT to be the cause or cure.

### Known limitation: cone/base corner numerical instability
A persistent solver instability originates at the cone/base corner (x≈0.4166m, r≈0.1494m), causing solver failure (`sigFpe` / iteration-limit abort) at a consistent simulation time of ~1.62–1.65×10⁻⁸ s across multiple mesh/timestep configurations. Root cause not established. Full investigation, evidence, and rejected/attempted fixes are documented in CASE_SPECIFICATION.md. **This must be revisited if downstream validation results (Fay-Riddell, Billig) show unphysical behavior without other explanation.**

### Status
Baseline mesh (short-domain, Case 4) is topologically verified and considered the primary Project 01 baseline going forward. Solver-viability testing revealed the above instability, which is deprioritized per explicit decision, not fixed. Next planned step: grid convergence study (verification) and validation against Fay-Riddell/Billig correlations, proceeding with awareness of this open item.

---

## UPDATE — 2026-09-08: Numerical Sensitivity Testing and Fan-Point Mesh Experiment

### Summary
Following the corner-instability investigation documented above, this session tested numerical-dissipation sensitivity (PIMPLE outer-correction count, convection scheme) and a targeted mesh-topology fix (Gmsh BoundaryLayer FanPointsList) against the known cone/base corner instability. The mesh fix produced a dramatic (~28x) improvement in survival time but did NOT eliminate the instability — it relocated to an adjacent, topologically clean region. The instability is now understood to be more likely a genuine local flow-physics challenge at the sharp convex corner than a pure mesh-topology artifact, though mesh sensitivity remains a contributing factor.

### Key results
1. **PIMPLE outer-correction sensitivity (nOuterCorrectors 1→10):** delayed failure by ~9.6% (1.624e-8s → 1.781e-8s). Confirms vanilla `shockFluid` does not deeply exploit outer correction for stability (consistent with literature comparison to shockFluidX, which explicitly adds this capability).
2. **Global scheme sensitivity (vanAlbada→Minmod):** delayed failure by ~18.4% (1.624e-8s → 1.924e-8s). Same failure mechanism, same corner location. Confirms scheme dissipation has a real but insufficient effect.
3. **Fan-point mesh fix (adding `p_base_outer` to `Field[1].FanPointsList`):** Deep topology diagnosis of Cells 4632-4637 revealed the true anomaly was not aspect ratio (present throughout wall_base, including in perfectly healthy cells) but a severely skewed (skew≈0.65) face connecting the wall_base BL stack directly to the wall_cone BL stack's terminal cell. Adding the corner point to FanPointsList (Gmsh's documented mechanism for exactly this two-BL-edges-meeting-at-a-point scenario) resolved this specific skew (→0.0000-0.0031, matching healthy reference cells) with zero cost to aspect ratio, non-orthogonality, or other settings.
4. **Result: failure delayed from 1.624e-8s to ~4.32-4.52e-7s — approximately 28x longer survival**, with the original epicenter cells (4632-4637) now stable and well-behaved. However, the instability re-emerged one radial layer inward (Cells 19290-19294), in cells confirmed via the same topology diagnostic to be PERFECTLY clean (0.00° non-orthogonality, 0.0000 skewness) — ruling out a topological mesh explanation at the new location.
5. **Time-history tracking of the corner region confirms the flow was still actively, monotonically evolving (not quasi-steady) at the time of failure** — meaning even the fan-point mesh's current state is not yet suitable for grid-convergence work without further investigation.

### Revised understanding
The corner instability is most likely driven by genuine, severe local flow physics (a sharp convex expansion at the cone/base corner) that is difficult for this solver/scheme/mesh combination to resolve robustly — mesh topology fixes can significantly delay and relocate the failure but have not eliminated it. This is a stronger, more specific characterization than the previous "known limitation, cause unresolved" framing.

### Status
Baseline mesh remains the adopted short-domain Case 4 configuration (fan-point fix NOT yet promoted to the primary baseline — remains an experimental variant, `openfoam/baseline_case_shortdomain_fanpoint_test/`). Grid convergence work is deferred pending either (a) further investigation into achieving genuine quasi-steady behavior before failure, or (b) an explicit decision to proceed with validation using only quantities that stabilize well before the failure point, if any do.

---

## UPDATE — 2026-09-08 (cont.): Corner/Base Investigation Concluded — Modeling-Domain Limitation

### Summary
The cone/base corner instability investigation is concluded. Root cause: a portion of the base/wake region reaches local Knudsen numbers far beyond continuum-breakdown thresholds (Kn up to ~400, vs. the standard 0.05-0.1 breakdown criterion), meaning the perfect-gas continuum Navier-Stokes model is not physically valid there -- independent of mesh, scheme, or timestep, all of which were tested and ruled out as the primary cause. This is a spatially LOCALIZED modeling-domain limitation, not a defect invalidating the full solution. Full technical detail in CASE_SPECIFICATION.md.

### Decision: Path 3A adopted
Project 01 proceeds within its original locked-in scope (perfect-gas, continuum CFD). The base/wake rarefaction is documented as an out-of-scope region for this project. A rarefied/DSMC or hybrid treatment, which would be needed to resolve it, is noted as candidate scope for a future, separate project -- not a Project 01 modification.

### Revised QoIs
Primary (proceeding to validation): stagnation-point heat flux, shock stand-off distance, forebody pressure distribution.
Deferred: total drag, pending quantification of base-pressure contribution to axial force.

### Status
Proceeding to temporal-independence study of forebody QoIs (using existing fan-point run data) before beginning grid convergence.
