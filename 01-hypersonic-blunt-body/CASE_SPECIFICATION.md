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

**Shock-layer refinement zone (Option C — broad bounded region, no a priori shock-location guess):**

A moderately fine, broadly bounded refinement zone is applied around the body
(Gmsh `Box` field) rather than either (a) a narrow box hand-placed at a
predicted shock location, or (b) solution-adaptive mesh refinement (AMR).
This avoids the risk of under-resolving a shock that lands outside a
narrowly-guessed region, while remaining simpler and fully reproducible
from the `.geo` file alone, without a solve-refine iteration loop. AMR is
left as a possible future upgrade if mesh-convergence verification shows
this is insufficient.

All parameters are expressed in terms of R_n, L, R_b so this structure
generalizes across the nose-radius sweep without modification.

| Parameter | Formula | Baseline value |
|---|---|---|
| Zone upstream edge | −3 × R_n | −0.15 m |
| Zone downstream edge | L (full body length) | 0.4166 m |
| Zone radial extent | R_b + 4 × R_n | 0.35 m |
| Refined cell size | R_n / 12 | ≈4.17 mm |
| Transition (blend) thickness | 2 × R_n | 0.10 m |

Verified (baseline case): 39,067 nodes, 38,893 triangles, 19,237 quads
(quads concentrated in the boundary-layer region). Visual and quantitative
checks confirm the boundary-layer field and shock-layer Box field coexist
correctly without disrupting each other.

## Solver

- `rhoCentralFoam` (OpenFOAM 11), steady-state operation via appropriate
  pseudo-transient/local-time-stepping approach or steady-state analog
  (exact time-stepping/convergence strategy to be finalized at solver-setup
  stage and documented there).

## Status

Baseline and sweep geometry parameters confirmed (sweep upper bound corrected
from 0.20 m to 0.14 m after tangency-validity check). Geometry generation not
yet started.

---

## UPDATE — 2026-09-06: Corner Defect Resolution, Domain Restructuring, and Corner Instability Investigation

### 1. Cone/Base Corner Topology Defect — RESOLVED

**Defect:** `gmshToFoam`-converted baseline wedge mesh showed 2 open cells and 1 misoriented face at the cone/base corner (x≈0.416, r≈0.149–0.150).

**Root-cause investigation:** Four hypotheses tested and rejected (shock-layer refinement zone boundary; missing BL fan-point treatment; local point sizing; sharp geometric corner itself via a rigorously-derived and verified 3mm fillet — all showed zero or negative effect on the defect).

**Confirmed cause (2×2 controlled experiment):**

| Case | Corner geometry | BL field `EdgesList` | Open cells | Misoriented faces |
|---|---|---|---|---|
| 1 | Sharp | {2,3} (original) | 2 | 1 |
| 2 | Filleted (3mm) | {2,3} (original) | 2 | 1 |
| 3 | Filleted (3mm) | {2,3,9,4} (extended) | 0 | 0 |
| 4 | Sharp | {2,3,4} (extended) | 0 | 0 |

**Conclusion:** The BoundaryLayer field's `EdgesList` topology transition — not wall geometry — caused the defect. **Case 4 (sharp corner + extended EdgesList) adopted as primary baseline.** Fillet confirmed unnecessary. All four cases and derivation files preserved (`mesh/case1_...`, `mesh/case3_...`, `mesh/archive_...`).

**New observation on Case 4 (unrelated, monitored):** One highly skew face (skewness 4.00243) found in the wake region (x≈0.551, r≈0.005), ~14cm aft of the base. Confirmed unrelated to the corner fix. Resolved incidentally by the domain-shortening change below (excluded from the shortened mesh; post-shortening max skewness = 3.75573, OK).

### 2. Domain Restructuring — Outlet Moved to L + 0.3×R_n

**Motivation:** Solver-viability testing on the clean Case 4 baseline (full 3×L wake domain) crashed after only 21 timesteps — fewer than on the pre-fix defective mesh (45 steps) — with `sigFpe` traced to genuine negative absolute pressure (p_min = -13.38 Pa) near the base corner/wake, distinct from the resolved topology defect.

**Decision:** Shorten downstream domain to remove most of the wake region, per explicit approval, accepting the trade-off against the original wake-clearance rationale (outlet-BC-artifact risk), since the wake region was already a documented low-fidelity limitation.

**Implementation (3 attempts):**
- Attempt 1 (outlet exactly at x=L): FAILED — `wall_base` and `outlet` became exactly collinear at `p_base_axis`, a degenerate 180° vertex Gmsh's mesher cannot resolve.
- Attempt 2 (outlet at x=L+0.00005): FAILED, worse — near-parallel curves over a finite arc length produced genuine spurious self-intersections (264 errors).
- **Attempt 3 (ADOPTED): short horizontal axis stub (`x_outlet = L + 0.3*R_n ≈ 0.4316m`) + genuine right-angle vertical outlet.** Restored original 8-curve topology (`Curve Loop(1) = {1,2,3,4,5,6,7,8}`), just with curve 5 (axis segment) shortened from ~0.83m to ~1.5cm. Verified clean in both `baseline_mesh_full.geo` and `baseline_mesh_wedge.geo`.

**Resulting short-domain mesh (post wedge-retype) verification:**

| Metric | Long-domain (3×L) | Short-domain (L+0.3Rn) |
|---|---|---|
| Cells | 54,335 | 48,851 |
| Max aspect ratio | 761.778 OK | 761.778 OK (identical) |
| Max skewness | 4.00243 (flagged) | 3.75573 (OK — wake skew face excluded) |
| Non-orthogonality | max 89.99°, avg 9.13° | max 89.92°, avg 8.58° (comparable) |
| Wedge planarity | benign FP warning (~2.5e-8m) | benign FP warning (~2.8e-8m, same class) |
| Failed checks | 2 | 1 |

**Short-domain mesh is the new primary baseline**, strictly cleaner than the prior long-domain mesh. Located: `mesh/baseline_mesh_wedge.geo`/`.msh`, case: `openfoam/baseline_case_shortdomain/`.

**Note:** A transient pre-retype aspect-ratio flag (1045.12 vs the true 761.778) was traced to an artifact of computing `checkMesh`'s aspect-ratio metric before wedge patches were retyped — resolved automatically upon retyping; not a real geometric difference. Mechanism not fully explained; noted as a minor unresolved curiosity, not consequential.

### 3. Cone/Base Corner Numerical Instability — UNRESOLVED, ACCEPTED AS KNOWN LIMITATION

**Symptom:** Solver (`rhoCentralFoam`) fails via `sigFpe` in `hePsiThermo::calculate()` (or, in one variant, a Ux-solver iteration-limit abort feeding into the same energy/temperature failure) at a highly consistent simulation time across configurations:

| Configuration | Failure time |
|---|---|
| Long-domain (3×L), maxCo=0.5 | 1.62391e-08 s |
| Short-domain (L+0.3Rn), maxCo=0.5 | 1.62449e-08 s |
| Short-domain, maxCo=0.1 (5x tighter) | ~1.6524e-08 s |

**Controlled experiments performed (all diagnostic-only, no fixes applied):**

1. **Domain-shortening test:** Removing the downstream wake did not change failure time or location. **Rules out wake/outlet as the cause.**
2. **Timestep-tightening test (maxCo 0.5→0.1, maxDeltaT 1e-6→1e-8):** Failure still occurred at essentially the same simulation time despite 5x tighter Courant control (actual Courant number stayed well under limit, ~0.10–0.12, throughout). **Rules out pure timestep under-resolution as the primary cause.**
3. **Spatial diagnosis (direct field inspection at last pre-crash timestep, short-domain run):** Localized collapse identified precisely at Cells 4632–4635 (x≈0.4166–0.4166, r≈0.1468–0.1494) — the cone/base corner. p drops to 131.5 Pa (vs ~5900–13800 Pa upstream), rho to 3.97e-4 kg/m³ (vs ~0.025–0.038 upstream), T elevated to 627–1270K in this cluster, Ux collapses from ~2100 m/s to 216–770 m/s across a ~17μm axial span. T remains near-freestream elsewhere; the collapse is NOT a distributed/global phenomenon.
4. **Local neighborhood diagnostic (BFS depth-2 around Cells 4632/4633):** Confirmed collapse spans a multi-cell cluster (not a single-cell pathology), forming a coherent low-p/rho patch precisely at the point where fine cone-wall tangential spacing meets the newly-BL-controlled base wall (per the Case 4 fix). Every cell in this cluster shows bbox aspect ratio 400–870, vs 12–20 for normal BL cells elsewhere (including a matched reference cell far from the corner). **Classification: (C) — a broader boundary-layer/corner interaction, not an isolated single-cell numerical artifact.** Whether the pressure collapse reflects genuine sharp-corner expansion physics (amplified by extreme cell geometry) versus a pure numerical artifact of that geometry remains **unresolved and unproven** — flagged explicitly, not assumed either way.
5. **Isolated sub-domain reconstruction attempt (topoSet/subsetMesh extraction of the corner region + parent-snapshot field mapping via nearest-neighbor):** Sub-mesh extraction and verification succeeded cleanly (topology, aspect ratio, patch splitting all confirmed valid). However, the reconstructed initial condition (independent per-cell nearest-neighbor mapping of T/p/rho/U with no cross-field consistency enforcement) failed almost immediately (2 timesteps, t=1.19e-9s) — 14x earlier than the parent runs. **Inconclusive**: this reflects a limitation of the naive field-reconstruction methodology (no access to the parent's actual `e` field, no EOS-consistency enforcement across independently-mapped fields), not new evidence about the corner's actual physics. OpenFOAM's `mapFields` utility was identified as the more appropriate tool for this if revisited. Sub-case preserved at `openfoam/corner_isolation_test/` as a documented negative result.
6. **Attempted mitigation — separate per-wall BoundaryLayer growth parameters:** Investigated whether Gmsh supports per-edge `hwall_n`/`ratio` within one BoundaryLayer field (it does not — confirmed via documentation search across CFD-Online and Gmsh's GitLab issue tracker: these are field-level properties applied uniformly across the field's entire `EdgesList`). Attempted workaround: two separate `BoundaryLayer` fields (fine on nose/cone, coarse on base wall) combined via a `Min` field, per a documented pattern from Gmsh's own issue tracker. **Result: NOT VIABLE in Gmsh 4.12.1.** The `Min` combination correctly merges scalar mesh-sizing values but silently discards structured BL quad-generation behavior entirely (verified: 0 quad elements in the resulting mesh despite `Quads=1` on both constituent fields; falls back to plain unstructured triangulation). An initial attempt using `BoundaryLayer Field = <Min field index>` segfaulted; removing that directive avoided the crash but also eliminated the intended structured-BL behavior. This mechanism is therefore not usable for achieving differentiated per-wall-segment BL resolution in this Gmsh version.

**DECISION (explicit, by user):** The corner instability is **accepted as a known, unresolved limitation** of the current mesh methodology. No further active debugging at this time. **Trigger for revisiting:** if grid convergence, validation (Fay-Riddell/Billig), or any downstream quantity-of-interest shows unphysical behavior without another clear explanation, this corner instability is the prime suspect and should be revisited first — including reconsidering the previously-shelved 3mm fillet (tested only against the now-resolved topology defect, never against this distinct instability) and/or a proper `mapFields`-based sub-domain reconstruction.

**Files preserved from this investigation:** `openfoam/baseline_case_shortdomain_tighttimestep/` (timestep-sensitivity test case), `openfoam/corner_isolation_test/` (sub-domain isolation attempt, includes `topoSetDict`, `createPatchDict`), `scripts/mesh/diagnose_crash_location_generalized_v2.py`, `scripts/mesh/local_neighborhood_diagnostic.py`, `scripts/mesh/build_isolated_subcase_fields.py`, `scripts/mesh/compare_corner_cell_topology.py`, `scripts/mesh/inspect_cell_direct.py`, `mesh/test_multi_bl_field.geo` (Gmsh Min-combination viability test, negative result).
