# Wedge Mesh Generation Methodology (Validated)

This document records the validated procedure for converting a 2D (x, r)
axisymmetric Gmsh surface mesh into an OpenFOAM axisymmetric wedge mesh
compatible with `rhoCentralFoam`. Established via Stage 4a viability testing
on a simplified test geometry (`wedge_rotate_v2.geo`), confirmed to produce
a solver-runnable mesh (74+ successful `rhoCentralFoam` timesteps at full
M=7 baseline conditions).

## Summary of the validated pipeline

1. Build the 2D (x, r) surface in Gmsh as usual (points, splines/lines,
   curve loop, plane surface).
2. `Recombine Surface{1};` on the 2D surface BEFORE extrusion, forcing a
   quad-dominant (via Gmsh's Blossom algorithm) surface mesh.
3. Rotate the surface to `theta = -half_angle` FIRST (`Rotate {...}
   { Surface{1}; }`), then rotationally extrude by the FULL angle
   (`Extrude {...} { Surface{1}; Layers{1}; Recombine; }`), sweeping from
   `-half_angle` to `+half_angle`. This produces a symmetric wedge in one
   step and preserves exact point correspondence between front/back faces
   (critical - this is what actually avoids OpenFOAM wedge-planarity
   failures, unlike the alternative front/back-patch-based
   `extrudeMesh`-utility approach, which repeatedly crashed - see
   "Rejected approach" below).
4. Tag `Physical Surface` for all boundaries EXCEPT curves that lie
   entirely on the rotation axis (r=0 at both endpoints) - such curves
   produce NO swept surface under rotation (a curve on the axis is
   invariant under rotation about that same axis) and will NOT have a
   corresponding `out[]` index. Attempting to tag them causes
   "Uninitialized variable" Gmsh errors.
5. Export via `-format msh2` (MSH 2.2 ASCII) - `gmshToFoam` in OpenFOAM 11
   does not reliably parse Gmsh's newer default MSH 4.1 format.
6. `gmshToFoam <file>.msh` to convert.
7. The un-tagged axis-region faces land in `defaultFaces`. Retype
   `frontWedge`/`backWedge` from generic `patch` to `wedge` type via DIRECT
   TEXT EDITING of `constant/polyMesh/boundary` (NOT `createPatch`, which
   was found to only reorder faces of an existing same-named patch without
   changing its type).
8. Run `checkMesh`. Expect (and accept) two categories of benign residual
   warnings:
   - Wedge patch non-planarity on the order of 1e-8 to 1e-15 metres
     (floating-point precision artifact of the strict internal check -
     matches a documented, known-benign OpenFOAM issue, OpenFOAM bug
     tracker #0001486).
   - A handful of degenerate/near-zero-area faces immediately adjacent to
     the rotation axis (in `defaultFaces`), typically <0.1% of total faces,
     concentrated at geometric singularities (e.g. a nose tip touching the
     axis). Multiple remediation attempts (`empty` patch retyping,
     `collapseEdges` with broad and targeted face sets, `Coherence`) did
     NOT resolve this cleanly and in two cases caused `collapseEdges` to
     crash with an internal `edgeCollapser::syncCollapse` error - this
     appears to be a genuine limitation/edge case in this OpenFOAM 11
     build's mesh-filtering tools for this specific topology, not a fixable
     configuration issue.
9. CONFIRMED: despite the above residual warnings, `rhoCentralFoam` (with
   `momentumTransport` set to laminar, standard PIMPLE/diagonal solver
   settings) runs successfully through many timesteps without any
   mesh-related fatal error. The defects do not block the solver.

## Rejected approach: OpenFOAM's `extrudeMesh` utility

An alternative approach - building a thin flat "slab" mesh in Gmsh (straight,
non-rotating extrude), converting via `gmshToFoam`, then using OpenFOAM's
own `extrudeMesh` utility (`extrudeModel wedge`, `constructFrom patch`) to
perform the actual wedge revolution - was attempted first and REJECTED
after repeated, unresolved crashes (`FOAM FATAL ERROR` / floating-point
exception inside `wedgePolyPatch::calcGeometry`, "not marked for collapse"
in `edgeCollapser::syncCollapse`). Root cause was never fully established;
suspected face-winding/orientation inconsistency introduced during
`gmshToFoam`'s reconstruction of patch face connectivity from the flat
slab's volume-cell topology, persisting even after eliminating tetrahedral
cells (switching to a clean hex-only slab via `Recombine Surface`) and
after using Gmsh's `Coherence` command. This approach is NOT recommended
for this project; the direct Gmsh-native rotate-then-extrude method
(above) should be used instead.

## Outstanding minor items (documented as accepted limitations)

- The near-axis degenerate face defect (defaultFaces, <0.1% of faces) is
  NOT fully resolved, only confirmed non-blocking for the solver. This
  will be re-examined if it causes any issue during the real baseline
  mesh solve or during post-processing near the stagnation point
  specifically (where peak heat flux is measured - a region of particular
  importance for this project's research question).
- Solver dictionary completeness (fvSchemes/fvSolution keyword
  requirements: `laplacianSchemes`, `PIMPLE`, `XFinal` variants for all
  solved fields, correct preconditioner-vs-matrix-symmetry pairing) was
  established through iterative trial against OpenFOAM 11's actual error
  messages, not from a single authoritative reference. The dictionaries in
  this test case should be treated as a validated STARTING TEMPLATE, not
  yet a fully tuned production solver configuration (e.g. deltaT was a
  rough placeholder with no CFL-based adjustTimeStep control, which is
  almost certainly why the test run eventually crashed with a
  floating-point exception in thermodynamic property calculation after
  ~74 steps - this is a numerics/timestep issue, unrelated to mesh
  topology, and will be properly addressed with adjustTimeStep/maxCo when
  building the real production solver setup).
