# Wall_cone/Wall_base Corner Mesh Defect — Investigation Log

## Symptom (stable across all attempts below)

`checkMesh` on the axisymmetric wedge mesh consistently reports:
- 2 "open cells" (max cell openness up to 1.0) at approximately (x≈0.416, r≈0.149)
- 1 face with incorrect orientation
- Located at the wall_cone / wall_base junction (aft end of the body,
  x≈L=0.41662, r≈R_b=0.15)

When solved with `rhoCentralFoam` (full M=7 baseline conditions), this
defect produces a confirmed non-physical temperature (T=196.4K, below
freestream T∞=226.5K, at the exact defective cell) after ~45 timesteps
under stable, CFL-controlled conditions, followed by a `sigFpe` crash in
`hePsiThermo::calculate()`. Applying `limitTemperature` (fvConstraint,
T bounded to [150,2500] K) did NOT prevent the crash - it still occurred
at a similar step count, indicating the numerical corruption originates
in the solved field (e or rho) before T is even computed, not something
a T-bound alone can catch.

## Four hypotheses tested, all REJECTED (no measurable effect on the defect)

### 1. Shock-layer refinement zone boundary coincident with the corner
**Hypothesis:** `Field[2].XMax = L` placed the Box field's fine/coarse
mesh-density transition exactly at the sharp geometric corner, causing
locally ambiguous cell formation.
**Fix attempted:** `Field[2].XMax = L + 0.1*R_n` (small margin past base).
**Result:** Identical defect (2 open cells, 1 misoriented face, same
non-orthogonality count). No measurable change.

### 2. Missing fan-point treatment at the corner
**Hypothesis:** The `BoundaryLayer` field's structured quad growth needs
explicit `FanPointsList` treatment at sharp corners (this worked for the
nose-tip/axis corner earlier in this project).
**Fix attempted:** `Field[1].FanPointsList = {p_nose_tip, p_base_outer}`.
**Result:** Identical defect. Zero effect - not even a change in mesh
node/element count, suggesting this setting did not influence this
region of the mesh at all.

### 3. Local point sizing (sliver-face relief)
**Hypothesis:** A genuine ~150-200x undersized sliver face was found at
this location (face 38437, area 3.3e-08 vs neighbouring faces ~5e-6) -
possibly from two near-coincident points.
**Fix attempted:** `p_base_outer` given `lc_wall * 3` (locally coarser
point spacing).
**Result:** Defect persisted (2 open cells, max openness now exactly 1.0
rather than 0.999999 - arguably marginally worse). Fix also had a large,
unintended side effect: wall_base and wall_cone face counts dropped
substantially (unintended broad coarsening), for no benefit. Reverted.

### 4. Sharp geometric corner itself (3mm fillet)
**Hypothesis:** A genuine C0 discontinuity in wall-normal direction at
the sharp cone/base corner is structurally incompatible with clean
per-face extrusion under wedge rotation.
**Fix attempted:** Explicit 3mm fillet (r_fillet = 0.003 m), derived
analytically (tangent-length/centre-offset trigonometry), verified to
machine precision (tangency dot products ~1e-18, radii exact to 0.003000m
at both tangent points), and visually confirmed smooth/continuous at
fine mesh resolution. Full engineering justification documented
(negligible relative to R_n, R_b; located far from primary QoIs -
stagnation heating, cone-flank Cp/heat-flux, shock stand-off).
**Result:** Identical defect signature (2 open cells, max openness 1,
1 misoriented face). This is the strongest negative result of the four:
a real, verified geometry change with zero measurable effect on the
mesh-topology defect strongly suggests the defect is NOT caused by wall
geometry at all.

## Leading hypothesis (NOT YET CONFIRMED)

`Field[1].EdgesList = {2, 3}` (the BoundaryLayer field's structured-quad
region) has been fixed at "nose curve + cone curve only" through ALL
FOUR fix attempts above, including the fillet (where curve 3 now
terminates earlier, at the fillet's tangent point, but EdgesList was
never extended to include the new fillet curve 9 or the base curve 4).

This means the structured-quad-to-unstructured-triangle mesh TOPOLOGY
transition has been pinned at the exact same physical wall location
(end of original curve 3) throughout every experiment, regardless of
what wall geometry sits downstream of it. This is consistent with the
defect's complete insensitivity to all four geometry/field changes
tested so far, all of which altered the WALL SHAPE but never the
BoundaryLayer field's own topology extent.

**This is a hypothesis, not a confirmed root cause.** It has not yet
been tested directly. The next step (per explicit instruction) is a
controlled topology experiment: extend `Field[1].EdgesList` to include
the fillet/base curves (or deliberately move this boundary) and observe
whether the defect moves or disappears, WITHOUT any further geometry
changes - isolating the topology-transition variable specifically.

## Preserved states

- `mesh/archive_sharpcorner_baseline_mesh_full.geo` - original sharp
  corner, 2D reference (retrieved from git commit aa9422f)
- `mesh/archive_filleted_baseline_mesh_full.geo` - current 3mm filleted
  2D reference
- `mesh/archive_filleted_baseline_mesh_wedge.geo` - current 3mm filleted
  wedge-extruded version (the one tested against checkMesh/solver above)

## Explicit non-action

Aggressive `T`/`e`/`rho` bounding has been deliberately NOT adopted as a
substitute fix. `limitTemperature` was tested as a DIAGNOSTIC SAFEGUARD
only (per earlier instruction) and confirmed NOT to prevent the crash -
this result itself is evidence the defect is a real numerical/topology
problem, not merely an extreme-but-otherwise-valid state that bounding
could reasonably contain.

## FINAL: Controlled 2×2 experiment matrix and causal conclusion

| Case | Corner geometry | BoundaryLayer EdgesList | Open cells | Misoriented faces | Result |
|---|---|---|---|---|---|
| 1 | Sharp | {2,3} (original) | 2 | 1 | DEFECTIVE (original baseline) |
| 2 | Filleted (3mm) | {2,3} (original) | 2 | 1 | DEFECTIVE (fillet alone: no effect) |
| 3 | Filleted (3mm) | {2,3,9,4} (extended) | 0 | 0 | CLEAN |
| 4 | Sharp | {2,3,4} (extended) | 0 | 0 | CLEAN |

**Causal conclusion:** The evidence strongly supports the BoundaryLayer
field's `EdgesList` topology transition as the cause of the observed
open-cell/misoriented-face defect, independent of wall geometry. Comparing
Case 1→2 (geometry changed, EdgesList unchanged: defect persists
identically) against Case 2→3 and Case 1→4 (EdgesList changed: defect
fully eliminated in both sharp and filleted geometries) isolates EdgesList
as the operative variable. The 3mm fillet is confirmed NOT REQUIRED to
resolve this specific defect.

**Adopted baseline: Case 4** (original sharp-corner geometry + extended
EdgesList {2,3,4}) - promoted to `mesh/baseline_mesh_full.geo` and
`mesh/baseline_mesh_wedge.geo`. This preserves the original, unmodified
body geometry (no engineering simplification needed) while resolving the
mesh defect through correct BoundaryLayer field configuration.

Case 3 (filleted + extended EdgesList) is preserved as a comparison/
sensitivity case: `mesh/case3_filleted_extendedEdges_wedge.geo/.msh`.
Case 1 (sharp + original EdgesList, the original defective configuration)
is preserved for the record: `mesh/case1_sharpcorner_originalEdges_wedge.geo`.

## New observation (Case 4 checkMesh, distinct from the resolved defect)

Case 4's `checkMesh` reports one remaining quality flag NOT present in
Case 3: `Max skewness = 4.00243, 1 highly skew faces detected`. Located
(via direct coordinate lookup, face ID 56102) at x≈0.551, r≈0.005 - in
the wake/axis-downstream region, ~14cm aft of the base (x=0.417), NOT
at the cone/base corner. Visually confirmed via
`scripts/mesh/plot_baseline_case4_check.py`
(`results/figures/baseline_case4_check.png`): the cone/base corner region
is clean and well-structured; the skew face sits in a dense band of thin
triangles along the centerline in the coarsening wake-transition zone,
plausibly from axis-line/surrounding-triangulation size mismatch - an
unrelated mechanism to the corner investigation.

Per explicit decision: this is recorded as a quantified, monitored
observation, NOT actioned with a geometry or mesh change at this stage.
It will be revisited only if solver-viability testing demonstrates
instability attributable to this location, or if it is shown to
materially affect a quantity of interest.

## STATUS: Mesh investigation phase CLOSED

Per explicit instruction, no further mesh-quality fixes will be pursued
at this stage. Next phase: controlled solver-viability testing of the
adopted Case 4 baseline (sharp-corner geometry, extended EdgesList),
inspecting numerical and physical behaviour before proceeding to
convergence/validation work.
