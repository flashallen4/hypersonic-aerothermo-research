#!/usr/bin/env python3
"""
Deep, read-only topology/geometry diagnostic of the cone/base corner cell
cluster (default seeds: 4632-4636) and immediate neighbors, plus off-corner
reference cells for comparison. Does NOT modify the mesh or run the solver.

Computes, per cell:
  - volume, bbox extent, bbox-based aspect ratio proxy
  - cell shape inferred from face count/type (hex/prism/pyramid/tet)
  - axial distance from wall_base plane (x - L), as a BL-layer-depth proxy
  - per-face: area, classification (boundary patch or internal+neighbor),
    and for internal faces: non-orthogonality angle and skewness computed
    directly from owner/neighbour centroids and face geometry.

Usage:
    python3 corner_topology_deep_diagnostic.py <case_dir>
"""

import re
import argparse
from pathlib import Path
import numpy as np


def find_balanced_block(text, start_idx):
    depth = 0
    i = start_idx
    while i < len(text):
        if text[i] == '(':
            depth += 1
        elif text[i] == ')':
            depth -= 1
            if depth == 0:
                return text[start_idx + 1:i], i + 1
        i += 1
    raise ValueError("Unbalanced parens")


def parse_foam_list_raw(path):
    text = Path(path).read_text()
    text = re.sub(r'//.*', '', text)
    text = re.sub(r'/\*.*?\*/', '', text, flags=re.DOTALL)
    m = re.search(r'\n\s*(\d+)\s*\n\s*\(', text)
    count = int(m.group(1))
    start = m.end()
    depth = 1
    i = start
    while depth > 0:
        if text[i] == '(':
            depth += 1
        elif text[i] == ')':
            depth -= 1
        i += 1
    return count, text[start:i - 1]


def parse_points(path):
    _, body = parse_foam_list_raw(path)
    pts = re.findall(r'\(([^)]+)\)', body)
    return np.array([[float(x) for x in p.split()] for p in pts])


def parse_faces(path):
    _, body = parse_foam_list_raw(path)
    faces = re.findall(r'\d+\(([^)]*)\)', body)
    return [[int(x) for x in f.split()] for f in faces]


def parse_labels(path):
    _, body = parse_foam_list_raw(path)
    return np.array([int(x) for x in body.split()])


def parse_boundary(path):
    text = Path(path).read_text()
    text = re.sub(r'/\*.*?\*/', '', text, flags=re.DOTALL)
    out = []
    for m in re.finditer(
        r'(\w+)\s*\{[^{}]*?type\s+(\w+);[^{}]*?nFaces\s+(\d+);[^{}]*?startFace\s+(\d+);[^{}]*?\}',
        text, re.S
    ):
        n, t, nf, sf = m.groups()
        out.append((n, t, int(nf), int(sf)))
    return out


def face_area_normal_centroid(pts):
    n = len(pts)
    centroid = pts.mean(axis=0)
    nv = np.zeros(3)
    for i in range(n):
        p1 = pts[i]
        p2 = pts[(i + 1) % n]
        nv += np.cross(p1 - centroid, p2 - centroid)
    area = 0.5 * np.linalg.norm(nv)
    unit_n = nv / (2 * area) if area > 0 else nv
    return area, unit_n, centroid


def classify_face(fid, n_internal, patches):
    if fid < n_internal:
        return "internal", None
    for name, ptype, nfaces, startface in patches:
        if startface <= fid < startface + nfaces:
            return "boundary", name
    return "boundary", "UNKNOWN"


def cell_shape_label(n_faces, face_vertcounts):
    quads = sum(1 for v in face_vertcounts if v == 4)
    tris = sum(1 for v in face_vertcounts if v == 3)
    if n_faces == 6 and quads == 6:
        return "hexahedron"
    if n_faces == 5 and quads == 3 and tris == 2:
        return "prism (triangular)"
    if n_faces == 5 and quads == 1 and tris == 4:
        return "pyramid"
    if n_faces == 4 and tris == 4:
        return "tetrahedron"
    return f"other (nfaces={n_faces}, quads={quads}, tris={tris})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--seed-cells", type=int, nargs="+", default=[4632, 4633, 4634, 4635, 4636])
    ap.add_argument("--bfs-depth", type=int, default=1)
    args = ap.parse_args()

    case = Path(args.case_dir)
    poly = case / "constant" / "polyMesh"

    # Geometry reference
    R_n, theta_c_deg, R_b = 0.05, 15.0, 0.15
    theta_c = np.radians(theta_c_deg)
    x_t = R_n * (1 - np.sin(theta_c))
    r_t = R_n * np.cos(theta_c)
    L = x_t + (R_b - r_t) / np.tan(theta_c)
    print(f"Reference: L (wall_base plane) = {L:.6f} m\n")

    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_labels(poly / "owner")
    neighbour = parse_labels(poly / "neighbour")
    patches = parse_boundary(poly / "boundary")
    n_internal = len(neighbour)

    cell_faces = {}
    for fid, o in enumerate(owner):
        cell_faces.setdefault(o, []).append(fid)
    for fid, n in enumerate(neighbour):
        cell_faces.setdefault(n, []).append(fid)

    def cell_centroid(cid):
        fids = cell_faces[cid]
        vids = set()
        for fid in fids:
            vids.update(faces[fid])
        pts = points[list(vids)]
        return pts.mean(axis=0)

    def cell_volume(cid):
        fids = cell_faces[cid]
        vol = 0.0
        for fid in fids:
            fpts = points[faces[fid]]
            area, unit_n, fcen = face_area_normal_centroid(fpts)
            sign = 1.0 if owner[fid] == cid else -1.0
            vol += sign * np.dot(fcen, unit_n) * area
        return vol / 3.0

    def report_cell(cid, label=""):
        fids = cell_faces[cid]
        vids = set()
        for fid in fids:
            vids.update(faces[fid])
        pts = points[list(vids)]
        bmin, bmax = pts.min(axis=0), pts.max(axis=0)
        ext = bmax - bmin
        nz = [e for e in ext if e > 1e-12]
        ar = max(nz) / min(nz) if nz else float('nan')
        cen = cell_centroid(cid)
        vol = cell_volume(cid)
        x, y, z = cen
        r = np.hypot(y, z)

        face_vertcounts = [len(faces[fid]) for fid in fids]
        shape = cell_shape_label(len(fids), face_vertcounts)

        n_internal_faces = sum(1 for fid in fids if fid < n_internal)

        print(f"\n{'='*70}")
        print(f"CELL {cid} {label}")
        print(f"{'='*70}")
        print(f"  centroid: x={x:.6f} r={r:.6f}  (x-L = {x-L:.6e} m)")
        print(f"  volume: {vol:.6e} m^3")
        print(f"  bbox extent: {ext}  (AR_bbox proxy = {ar:.1f})")
        print(f"  shape: {shape}  ({len(fids)} faces, neighbor_count={n_internal_faces})")
        print(f"  faces:")
        for fid in fids:
            fpts = points[faces[fid]]
            area, unit_n, fcen = face_area_normal_centroid(fpts)
            cls, pname = classify_face(fid, n_internal, patches)
            if cls == "internal":
                other = neighbour[fid] if owner[fid] == cid else owner[fid]
                Cp = cell_centroid(cid)
                Cn = cell_centroid(other)
                d = Cn - Cp
                dmag = np.linalg.norm(d)
                # non-orthogonality: angle between d and face normal
                cos_ang = np.dot(d, unit_n) / (dmag * 1.0) if dmag > 0 else 1.0
                cos_ang = np.clip(cos_ang, -1, 1)
                nonortho_deg = np.degrees(np.arccos(abs(cos_ang)))
                # skewness: distance from face centre to line P-N intersection with face plane, / |d|
                denom = np.dot(d, unit_n)
                if abs(denom) > 1e-30:
                    t = np.dot(fcen - Cp, unit_n) / denom
                    intersect = Cp + t * d
                    skew = np.linalg.norm(fcen - intersect) / dmag if dmag > 0 else 0.0
                else:
                    skew = float('nan')
                print(f"    face {fid}: area={area:.4e}  INTERNAL -> cell {other}  "
                      f"nonortho={nonortho_deg:.2f}deg  skew={skew:.4f}")
            else:
                print(f"    face {fid}: area={area:.4e}  BOUNDARY:{pname}")

    # BFS neighborhood
    def neighbors_of(cid):
        result = []
        for fid in cell_faces[cid]:
            if fid < n_internal:
                other = neighbour[fid] if owner[fid] == cid else owner[fid]
                result.append(other)
        return result

    visited = set(args.seed_cells)
    frontier = set(args.seed_cells)
    for _ in range(args.bfs_depth):
        new_frontier = set()
        for cid in frontier:
            for other in neighbors_of(cid):
                if other not in visited:
                    new_frontier.add(other)
        visited.update(new_frontier)
        frontier = new_frontier

    print(f"Neighborhood (BFS depth {args.bfs_depth} from seeds {args.seed_cells}): "
          f"{len(visited)} cells total\n")

    print("### SEED CELLS (detailed) ###")
    for cid in args.seed_cells:
        report_cell(cid, label="[SEED]")

    print("\n\n### FULL NEIGHBORHOOD SUMMARY TABLE ###")
    rows = []
    for cid in visited:
        cen = cell_centroid(cid)
        vol = cell_volume(cid)
        fids = cell_faces[cid]
        vids = set()
        for fid in fids:
            vids.update(faces[fid])
        pts = points[list(vids)]
        ext = pts.max(axis=0) - pts.min(axis=0)
        nz = [e for e in ext if e > 1e-12]
        ar = max(nz) / min(nz) if nz else float('nan')
        shape = cell_shape_label(len(fids), [len(faces[f]) for f in fids])
        x, y, z = cen
        r = np.hypot(y, z)
        wall_patches = set()
        for fid in fids:
            cls, pname = classify_face(fid, n_internal, patches)
            if cls == "boundary" and pname in ("wall_cone", "wall_base", "wall_nose"):
                wall_patches.add(pname)
        rows.append((cid, x, r, x - L, vol, ar, shape, ",".join(wall_patches) if wall_patches else "-"))
    rows.sort(key=lambda t: (t[1], t[2]))
    print(f"{'cell':>6} {'x':>10} {'r':>10} {'x-L':>12} {'vol':>12} {'AR':>8} {'shape':>20} {'wall'}")
    for cid, x, r, xL, vol, ar, shape, wall in rows:
        print(f"{cid:>6} {x:>10.6f} {r:>10.6f} {xL:>12.4e} {vol:>12.4e} {ar:>8.1f} {shape:>20} {wall}")

    # --- Off-corner reference cells for comparison ---
    print("\n\n### OFF-CORNER REFERENCE CELLS ###")
    wall_cone_patch = next((pp for pp in patches if pp[0] == "wall_cone"), None)
    wall_base_patch = next((pp for pp in patches if pp[0] == "wall_base"), None)

    def cell_owning_face(fid):
        return owner[fid]

    if wall_cone_patch:
        name, ptype, nfaces, startface = wall_cone_patch
        for frac, tag in [(0.5, "mid-cone"), (0.9, "near-corner-end-of-cone")]:
            fid = startface + int(nfaces * frac)
            cid = cell_owning_face(fid)
            report_cell(cid, label=f"[REFERENCE: wall_cone, {tag}]")

    if wall_base_patch:
        name, ptype, nfaces, startface = wall_base_patch
        for frac, tag in [(0.1, "wall_base near axis (far from corner)"),
                           (0.5, "wall_base mid-radius")]:
            fid = startface + int(nfaces * frac)
            cid = cell_owning_face(fid)
            report_cell(cid, label=f"[REFERENCE: {tag}]")


if __name__ == "__main__":
    main()
