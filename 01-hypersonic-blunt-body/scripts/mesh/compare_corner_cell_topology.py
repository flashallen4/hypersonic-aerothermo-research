#!/usr/bin/env python3
"""
Diagnostic (read-only): locate and compare the cone/base corner cell's
topology and geometry between two OpenFOAM cases (e.g. short-domain vs
long-domain baseline). Parses raw polyMesh ASCII files directly.
Does not modify anything.

Usage:
    python3 compare_corner_cell_topology.py <case_dir> [--x0 X --r0 R]

Run once per case directory; compare the two printed reports manually.
"""

import re
import sys
import argparse
from pathlib import Path
import numpy as np


def extract_block(text: str) -> str:
    idx = text.index('}')
    text = text[idx + 1:]
    start = text.index('(')
    depth = 0
    for i in range(start, len(text)):
        c = text[i]
        if c == '(':
            depth += 1
        elif c == ')':
            depth -= 1
            if depth == 0:
                return text[start + 1:i]
    raise ValueError("Unbalanced parentheses")


def read_points(path):
    text = Path(path).read_text()
    block = extract_block(text)
    entries = re.findall(r'\(([^()]+)\)', block)
    return np.array([[float(v) for v in e.split()] for e in entries])


def read_labels(path):
    text = Path(path).read_text()
    block = extract_block(text)
    nums = re.findall(r'-?\d+', block)
    return np.array([int(n) for n in nums])


def read_faces(path):
    text = Path(path).read_text()
    block = extract_block(text)
    faces = []
    for m in re.finditer(r'(\d+)\(([^()]+)\)', block):
        verts = [int(v) for v in m.group(2).split()]
        faces.append(verts)
    return faces


def read_boundary(path):
    text = Path(path).read_text()
    block = extract_block(text)
    patches = []
    for m in re.finditer(
        r'(\w+)\s*\{[^{}]*?type\s+(\w+);[^{}]*?nFaces\s+(\d+);[^{}]*?startFace\s+(\d+);[^{}]*?\}',
        block, re.S
    ):
        name, ptype, nfaces, startface = m.groups()
        patches.append((name, ptype, int(nfaces), int(startface)))
    return patches


def classify_face(face_id, n_internal, patches):
    if face_id < n_internal:
        return "internal"
    for name, ptype, nfaces, startface in patches:
        if startface <= face_id < startface + nfaces:
            return f"boundary:{name}"
    return "boundary:UNKNOWN"


def face_area_normal(pts):
    """Newell's method centroid/area for a planar-ish polygon."""
    n = len(pts)
    centroid = pts.mean(axis=0)
    normal = np.zeros(3)
    for i in range(n):
        p1 = pts[i]
        p2 = pts[(i + 1) % n]
        normal += np.cross(p1 - centroid, p2 - centroid)
    area = 0.5 * np.linalg.norm(normal)
    return area


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--x0", type=float, default=None, help="target x (default: cone/base corner, computed)")
    ap.add_argument("--r0", type=float, default=None, help="target r (default: R_b)")
    args = ap.parse_args()

    # Geometry reference (Project 01 locked baseline params)
    R_n = 0.05
    theta_c = np.radians(15.0)
    R_b = 0.15
    x_t = R_n * (1 - np.sin(theta_c))
    r_t = R_n * np.cos(theta_c)
    L = x_t + (R_b - r_t) / np.tan(theta_c)

    x0 = args.x0 if args.x0 is not None else L
    r0 = args.r0 if args.r0 is not None else R_b

    case = Path(args.case_dir)
    pm = case / "constant" / "polyMesh"

    points = read_points(pm / "points")
    faces = read_faces(pm / "faces")
    owner = read_labels(pm / "owner")
    neighbour = read_labels(pm / "neighbour")
    patches = read_boundary(pm / "boundary")
    n_internal = len(neighbour)
    n_cells = max(owner.max(), neighbour.max()) + 1 if len(neighbour) else owner.max() + 1

    print(f"Case: {case}")
    print(f"Target corner: x0={x0:.6f}, r0={r0:.6f}")
    print(f"Loaded: {len(points)} points, {len(faces)} faces, {n_internal} internal faces, ~{n_cells} cells")

    # --- find candidate cells near the corner by checking cell centroids ---
    def faces_of_cell(cell_id):
        fids = []
        for fid, o in enumerate(owner):
            if o == cell_id:
                fids.append(fid)
        for fid, n in enumerate(neighbour):
            if n == cell_id:
                fids.append(fid)
        return fids

    def cell_geom(cell_id):
        fids = faces_of_cell(cell_id)
        vert_ids = set()
        for fid in fids:
            vert_ids.update(faces[fid])
        pts = points[list(vert_ids)]
        centroid = pts.mean(axis=0)
        bbox_min = pts.min(axis=0)
        bbox_max = pts.max(axis=0)
        extent = bbox_max - bbox_min
        return centroid, bbox_min, bbox_max, extent, fids

    # Coarse scan: build a quick x,r index of face-based points to find nearby cells
    # We scan owner/neighbour cell ids whose faces have a vertex near (x0,r0)
    tol = 0.002  # 2 mm search radius in x, r
    candidate_cells = set()
    for fid, verts in enumerate(faces):
        pts = points[verts]
        # cylindrical r
        r_vals = np.hypot(pts[:, 1], pts[:, 2])
        x_vals = pts[:, 0]
        if np.any((np.abs(x_vals - x0) < tol) & (np.abs(r_vals - r0) < tol)):
            candidate_cells.add(owner[fid])
            if fid < n_internal:
                candidate_cells.add(neighbour[fid])

    print(f"\nCandidate cells near corner (within {tol*1000:.1f} mm): {len(candidate_cells)}")

    results = []
    for cid in sorted(candidate_cells):
        centroid, bmin, bmax, extent, fids = cell_geom(cid)
        x, y, z = centroid
        r = np.hypot(y, z)
        results.append((cid, x, r, extent, bmin, bmax, fids))

    # Sort by distance to target corner
    results.sort(key=lambda t: (t[1] - x0) ** 2 + (t[2] - r0) ** 2)

    print(f"\nTop 5 closest candidate cells to corner (x0={x0:.5f}, r0={r0:.5f}):")
    for cid, x, r, extent, bmin, bmax, fids in results[:5]:
        dx, dy, dz = extent
        # approximate aspect ratio proxy: largest extent / smallest nonzero extent
        nz_extent = [e for e in extent if e > 1e-12]
        ar_proxy = max(nz_extent) / min(nz_extent) if nz_extent else float('nan')
        print(f"  Cell {cid}: centroid=(x={x:.6f}, r={r:.6f})  "
              f"extent=(dx={dx:.6e}, dy={dy:.6e}, dz={dz:.6e})  AR_proxy={ar_proxy:.2f}")

    # --- Detailed report on best-match cell ---
    if results:
        cid, x, r, extent, bmin, bmax, fids = results[0]
        print(f"\n{'='*70}")
        print(f"DETAILED REPORT: Cell {cid} (closest match to corner)")
        print(f"{'='*70}")
        print(f"Centroid: x={x:.6f}, r={r:.6f}")
        print(f"Bounding box min: {bmin}")
        print(f"Bounding box max: {bmax}")
        print(f"Extent (dx, dy, dz): {extent}")
        print(f"Number of faces: {len(fids)}")
        print(f"\nFace-by-face breakdown:")
        for fid in fids:
            verts = faces[fid]
            pts = points[verts]
            area = face_area_normal(pts)
            cls = classify_face(fid, n_internal, patches)
            fx = pts[:, 0].mean()
            fr = np.hypot(pts[:, 1], pts[:, 2]).mean()
            neighbor_cell = None
            if fid < n_internal:
                neighbor_cell = neighbour[fid] if owner[fid] == cid else owner[fid]
            print(f"  face {fid}: area={area:.6e}  x={fx:.6f} r={fr:.6f}  "
                  f"{cls:25s} neighbor_cell={neighbor_cell}")


if __name__ == "__main__":
    main()
