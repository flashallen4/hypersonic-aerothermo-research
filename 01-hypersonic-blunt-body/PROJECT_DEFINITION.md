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

---

## Session Update: Axis-Topology Investigation, Shock Stand-Off Rejection, and Final Temporal Assessment

### Axis/Stagnation-Streamline Mesh-Topology Investigation

A dedicated mesh-topology diagnostic (`scripts/mesh/axis_topology_diagnostic.py`, read-only,
existing mesh only) was performed to resolve contamination discovered in an initial shock-detection
attempt.

**Finding:** the wedge mesh represents r=0 via 58 true axis points (points with r <= 1e-12) and 112
axis-adjacent faces, all of type `frontWedge`/`backWedge` (never internal or wall patches). These
resolve to exactly 56 unique axis-adjacent cells, forming a clean, monotonic, non-duplicated sequence
in x from far upstream (x=-0.28244 m) to the stagnation cell (48880, x=-0.00011 m) directly adjacent
to `wall_nose`. Zero near-duplicate-x pairs were found among these cells (uniqueness check).

By contrast, the previously-used r<0.003 m radius-tolerance selection pulled in ~7x more cells in the
same x-range (134 vs 18 in x in [-0.015, 0.005]), including cells at multiple distinct radii near
nearly identical x — an artifact of the body's curvature near the nose, not a true 1D axis profile.
In the bin [0.000, 0.001) specifically (the exact region of the previously-detected "shock" at
0.03-0.14mm), there were 0 true axis cells but 28 tolerance-selected cells, confirming that region is
entirely off-axis-cell contamination.

### Shock Stand-Off: Prior Detections Formally Rejected

Both prior shock-detection attempts are formally rejected as invalid:
- The original "first axis cell where p > 2*p_inf" method (result: ~0.11 mm stand-off).
- The subsequent peak-|dp/dx|/peak-|drho/dx| method applied to the r<0.003-tolerance-selected axis
  (result: ~0.03-0.14 mm stand-off, varying by snapshot).

Both were built on a contaminated cell selection (multiple radial rays interleaved at similar x),
producing large spurious gradients that do not correspond to a real physical shock crossing. This is
confirmed directly via the mesh-topology investigation above, not merely suspected.

### True-Axis Field Profile (uncontaminated, LATE snapshot t=4.48921e-07 s)

Using the 56 true axis-adjacent cells: clean freestream (p=1197 Pa, M=7.000) from x=-0.282 to
x=-0.012 m; a weak, small pre-shock disturbance (x=-0.0118 to -0.0092 m, p dips slightly, M rises
slightly above 7); then a genuine, smooth, monotonic compression from x=-0.008 m to the stagnation
cell (p: 1251 -> 1327 -> 1384 -> 1507 -> 1573 -> 1649 -> 1703 -> 1747 -> 1776 Pa; M: 6.84 -> 6.64 ->
6.49 -> 6.21 -> 6.06 -> 5.90 -> 5.80 -> 5.73 -> 5.68). One cell (48881, x=-0.00048 m) shows clearly
non-physical values (p=0.846 Pa, T=0.32 K, M=65.05) -- a numerical artifact, not a flow feature.
Stagnation cell 48880 shows p=12586.10 Pa, T=560.79 K, M=4.420, exactly matching the independently
computed p_stag value from the temporal-independence QoI script at the same snapshot -- confirming
internal consistency of the stagnation-cell identification.

**Conclusion (conservative, as required):** the nose flow contains a progressively strengthening
compression region, but a distinct, steady bow-shock discontinuity has not yet developed in the
currently available solution. The onset location (~8 mm upstream) is broadly consistent with the
Billig-correlation estimate for this geometry (R_n=0.05 m, M=7: delta ~= 7.6 mm), which is
encouraging, but the compression itself remains smeared/gradual (p reaches only ~1780-2040 Pa,
roughly 1.5-1.7x freestream, prior to the corrupted cell, and only ~10.5x freestream at the
stagnation cell itself) rather than exhibiting the near-instantaneous jump expected of a
well-captured M=7 normal shock (theoretical p2/p1=57.0, rho2/rho1=5.444, T2/T1=10.469). No claim is
made about whether this compression would eventually steepen into a sharp discontinuity given more
simulation time; the available data does not establish that either way. No new shock-detection
algorithm has been implemented as a result of this finding -- it is documented as a negative/open
result.

### Final Temporal Assessment (Primary Fan-Point Mesh, Frozen Configuration)

A final continuation run was performed strictly under the frozen configuration (geometry, mesh,
fvSchemes, solver settings, BCs all unchanged; only diagnostic output frequency, writeInterval,
had been previously tightened to 2 for data-capture purposes). No new intervention (no mesh/scheme/
solver/fvConstraints change) was introduced during this assessment.

**Dataset:** 24 snapshots, t=1.0604e-08 s to t=4.48921e-07 s.

**Convergence criterion (pre-defined):** |Q(i)-Q(i-1)|/Q(i-1) < 1% for 3 consecutive snapshots, with
a non-reversing (monotonically decreasing or flat) trend.

**A. q_stag temporal independence:** YES, by the formal criterion. Last three consecutive
   relative changes: 0.403%, 0.301%, 0.295% (all <1%, monotonically decreasing).

**B. p_stag temporal independence:** YES, by the formal criterion. Last three consecutive
   relative changes: 0.433%, 0.388%, 0.384% (all <1%, monotonically decreasing).

**C. Base-region breakdown vs QoI convergence:** breakdown occurred AFTER QoI convergence, but only
   marginally -- the convergence window closes at t=4.4892e-07 s; the deterministic sigFpe failure
   (same mechanism as previously documented: hePsiThermo::calculate()) occurs at approximately
   t=4.496-4.501e-07 s, roughly 7-10e-9 s (~10-14 iterations) later. This is characterized as a
   marginal pass, not a robust, far-from-breakdown plateau.

**D. Final usable temporal window / convergence times:** full usable window t=1.0604e-08 s to
   4.48921e-07 s; q_stag and p_stag convergence window t=4.4610e-07 s to t=4.4892e-07 s (same three
   snapshots for both quantities); failure at ~4.496-4.501e-07 s (exact value varies slightly,
   ~0.001-0.01e-07 s, across restart attempts due to ASCII checkpoint round-trip precision --
   consistent with previously documented restart-precision sensitivity, not a new finding).

**E. Final QoI status:**
  - Forebody surface pressure (5 wall_cone stations): temporally independent, converged early
    (~t=6.6e-8 s onward), stable ~8650-9200 Pa band, demonstrably unaffected by later base-region
    rarefaction growth. Strongest, most robust QoI in the project to date.
  - Stagnation heat flux (q_stag): meets the formal criterion, but only in a narrow (~2.8e-8 s)
    window close to breakdown -- classified as a MARGINAL PASS.
  - Stagnation pressure (p_stag): same -- MARGINAL PASS.
  - Shock stand-off distance: NO DEFENSIBLE VALUE. Both prior results (0.11 mm; 0.03-0.14 mm)
    formally rejected as axis-contamination artifacts (see above). True-axis analysis shows a
    progressively strengthening, not-yet-steady compression region -- documented as a genuine
    negative/open result, not a detection-method failure.
  - Drag: still deferred (base-pressure contribution to total axial force remains unquantified).

**F. Grid-convergence readiness:** Project 01 may proceed to grid convergence for FOREBODY SURFACE
   PRESSURE ONLY. Stagnation quantities (q_stag, p_stag) remain formally unresolved for grid-study
   purposes -- their marginal, near-breakdown convergence window is not considered a sufficiently
   robust basis, since mesh refinement is expected to alter the timescale on which base-region
   breakdown occurs, potentially shifting or eliminating this narrow window. Shock stand-off is not
   ready for grid convergence in any form; no baseline value exists to refine.

**G. No new intervention:** confirmed via direct git diff against the last committed
   thermophysicalProperties (zero difference) and via explicit grep of controlDict before and after
   this run. No mesh/scheme/solver/fvConstraints modification was made during this assessment. The
   final run's failure to capture one additional snapshot (falling ~10-14 iterations short of the
   writeInterval=2 threshold) was accepted as the natural boundary of the frozen configuration's
   available data, not treated as grounds for a further intervention cycle.

### New diagnostic scripts this session
- `scripts/mesh/shock_standoff_detection.py` -- peak-|dp/dx| / peak-|drho/dx| method; produced the
  now-rejected 0.03-0.14mm results using r<0.003-tolerance axis selection. Preserved for reference
  and as documentation of a rejected approach; NOT to be used for shock-location claims.
- `scripts/mesh/axis_structure_diagnostic.py` -- full axis-profile and multi-peak candidate-feature
  diagnostic (still uses r-tolerance selection; superseded for axis-cell selection purposes by the
  topology diagnostic below, but retained for its general profile/peak-characterization logic).
- `scripts/mesh/axis_topology_diagnostic.py` -- the decisive tool: identifies true axis-adjacent
  cells via degenerate-wedge-face connectivity (not radius tolerance), confirms a unique monotonic
  56-cell stagnation-streamline sequence, and is the basis for the topology conclusions above.
