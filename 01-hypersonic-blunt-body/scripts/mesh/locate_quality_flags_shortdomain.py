#!/usr/bin/env python3
"""
Diagnostic (read-only): locate the cell(s) in highAspectRatioCells and the
faces in nonOrthoFaces, written by checkMesh into constant/polyMesh/sets/.
Parses raw OpenFOAM ASCII polyMesh files directly (points, faces, owner,
neighbour, boundary, sets) via regex/bracket parsing -- no OpenFOAM Python
bindings required. Does not modify anything.

Run from the case directory, e.g.:
    cd openfoam/baseline_case_shortdomain
    python3 ../../scripts/mesh/locate_quality_flags_shortdomain.py
"""

import re
import sys
from pathlib import Path
import numpy as np

POLYMESH = Path("constant/polyMesh")
SETS = POLYMESH / "sets"

# --- Known Project 01 geometry parameters, for spatial classification ---
R_n = 0.05
theta_c_deg = 15.0
R_b = 0.15
theta_c = np.radians(theta_c_deg)
x_t = R_n * (1 - np.sin(theta_c))
r_t = R_n * np.cos(theta_c)
L = x_t + (R_b - r_t) / np.tan(theta_c)
x_outlet = L + 0.3 * R_n

REGIONS = [
    ("nose region",            -0.01,        x_t),
    ("cone wall region",        x_t,          L),
    ("corner / short axis stub", L - 0.01,    x_outlet + 0.01),
    ("outlet plane",            x_outlet - 0.005, x_outlet + 0.02),
]


def extract_block(text: str) -> str:
    """Return contents of the first top-level (...) block after the FoamFile header."""
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
    raise ValueError("Unbalanced parentheses while parsing block")


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


def region_label(x, r):
    for label, xmin, xmax in REGIONS:
        if xmin <= x <= xmax:
            return label
    if x < -0.01:
        return "farfield/upstream"
    if r > R_b + 4 * R_n - 0.05:
        return "farfield/outer"
    return "shock-layer/domain interior"


def main():
    print(f"Geometry reference values: x_t={x_t:.5f} L={L:.5f} x_outlet={x_outlet:.5f} R_b={R_b}")
    print()

    points = read_points(POLYMESH / "points")
    faces = read_faces(POLYMESH / "faces")
    owner = read_labels(POLYMESH / "owner")
    neighbour = read_labels(POLYMESH / "neighbour")
    patches = read_boundary(POLYMESH / "boundary")
    n_internal = len(neighbour)

    print(f"Loaded: {len(points)} points, {len(faces)} faces, "
          f"{len(owner)} owners, {n_internal} internal faces")
    print()

    # --- Build cell -> face list (only as needed, on demand) ---
    def faces_of_cell(cell_id):
        fids = []
        for fid, o in enumerate(owner):
            if o == cell_id:
                fids.append(fid)
        for fid, n in enumerate(neighbour):
            if n == cell_id:
                fids.append(fid)
        return fids

    def cell_centroid_and_bbox(cell_id):
        fids = faces_of_cell(cell_id)
        vert_ids = set()
        for fid in fids:
            vert_ids.update(faces[fid])
        pts = points[list(vert_ids)]
        centroid = pts.mean(axis=0)
        bbox_min = pts.min(axis=0)
        bbox_max = pts.max(axis=0)
        return centroid, bbox_min, bbox_max

    # --- 1. High aspect ratio cell ---
    print("=" * 70)
    print("HIGH ASPECT RATIO CELL(S)")
    print("=" * 70)
    har_path = SETS / "highAspectRatioCells"
    if har_path.exists():
        har_cells = read_labels(har_path)
        for cid in har_cells:
            centroid, bmin, bmax = cell_centroid_and_bbox(int(cid))
            x, y, z = centroid
            r = np.hypot(y, z)
            print(f"Cell {cid}: centroid = ({x:.6f}, {y:.6f}, {z:.6f})  "
                  f"[x={x:.5f}, r={r:.6f}]")
            print(f"  bbox min = {bmin}")
            print(f"  bbox max = {bmax}")
            print(f"  region classification: {region_label(x, r)}")
    else:
        print(f"File not found: {har_path}")
    print()

    # --- 2. Non-orthogonal faces ---
    print("=" * 70)
    print("SEVERELY NON-ORTHOGONAL FACES (>70 deg)")
    print("=" * 70)
    no_path = SETS / "nonOrthoFaces"
    if no_path.exists():
        no_faces = read_labels(no_path)
        print(f"Total faces in set: {len(no_faces)}")

        xs, rs, classes, regions = [], [], [], []
        for fid in no_faces:
            fid = int(fid)
            vert_ids = faces[fid]
            pts = points[vert_ids]
            cx, cy, cz = pts.mean(axis=0)
            r = np.hypot(cy, cz)
            xs.append(cx)
            rs.append(r)
            classes.append(classify_face(fid, n_internal, patches))
            regions.append(region_label(cx, r))

        xs = np.array(xs)
        rs = np.array(rs)

        print(f"x range: [{xs.min():.5f}, {xs.max():.5f}]")
        print(f"r range: [{rs.min():.5f}, {rs.max():.5f}]")
        print()

        print("Region breakdown (count per classified region):")
        from collections import Counter
        for region, count in Counter(regions).most_common():
            print(f"  {region:35s} : {count}")
        print()

        print("Internal vs boundary breakdown:")
        for cls, count in Counter(classes).most_common():
            print(f"  {cls:35s} : {count}")
        print()

        print("x-histogram (0.01 m bins):")
        bins = np.arange(xs.min() - 0.005, xs.max() + 0.015, 0.01)
        hist, edges = np.histogram(xs, bins=bins)
        for i, h in enumerate(hist):
            if h > 0:
                print(f"  x in [{edges[i]:.4f}, {edges[i+1]:.4f}) : {h} faces")
        print()

        print("First 10 sample face coordinates (id, x, r, class, region):")
        for i in range(min(10, len(no_faces))):
            fid = int(no_faces[i])
            print(f"  face {fid}: x={xs[i]:.6f} r={rs[i]:.6f} "
                  f"{classes[i]:25s} {regions[i]}")
    else:
        print(f"File not found: {no_path}")


if __name__ == "__main__":
    main()
