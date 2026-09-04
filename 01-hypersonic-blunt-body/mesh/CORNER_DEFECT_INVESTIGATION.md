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
