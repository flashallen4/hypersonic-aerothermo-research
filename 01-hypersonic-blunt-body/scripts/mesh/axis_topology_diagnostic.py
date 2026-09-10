#!/usr/bin/env python3
"""
Axis/stagnation-streamline mesh-topology diagnostic for Project 01
(read-only, existing mesh only; no mesh modification, no solver run).

Purpose: determine how the wedge mesh actually represents r=0 BEFORE
defining any shock-extraction algorithm. Does NOT assume "nearest centroid
to r=0" is the same as a true axis-adjacent cell. Does NOT assume the
wedge-plane representation gives a unique r=0 cell column -- establishes
this from actual point/face/cell connectivity.

Method:
  - True axis points: mesh points with r = hypot(y,z) exactly (or to
    floating-point precision) zero.
  - True axis-adjacent faces: faces containing >=2 axis points (i.e. an
    edge of the face lies exactly on r=0). For a 1-cell-thick wedge with
    a collapsed axis edge, such faces are typically degenerate (near-zero
    area) and live in defaultFaces/patch-type boundary, per prior mesh
    documentation -- but this script verifies that directly rather than
    assuming it.
  - True axis-adjacent cells: owner/neighbour cells of those faces.
  - For comparison, radius-tolerance cells (r < tol) are also computed,
    to explicitly show any difference between the two selections.
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


def face_type_for(fid, patches, n_internal):
    if fid < n_internal:
        return "internal"
    for name, ptype, nf, sf in patches:
        if sf <= fid < sf + nf:
            return f"{name}({ptype})"
    return "UNKNOWN_BOUNDARY"


def parse_scalar_field(path):
    if not path.exists():
        return None
    text = path.read_text()
    text_nc = re.sub(r'//.*', '', text)
    m = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\n?\s*(\d+)\s*\n', text_nc)
    if m is None:
        return None
    n = int(m.group(1))
    open_idx = text_nc.index('(', m.end() - 1)
    body, _ = find_balanced_block(text_nc, open_idx)
    vals = np.array([float(x) for x in body.split()])
    assert len(vals) == n
    return vals


def parse_vector_field(path):
    if not path.exists():
        return None
    text = path.read_text()
    text_nc = re.sub(r'//.*', '', text)
    m = re.search(r'internalField\s+nonuniform\s+List<vector>\s*\n?\s*(\d+)\s*\n', text_nc)
    if m is None:
        return None
    n = int(m.group(1))
    open_idx = text_nc.index('(', m.end() - 1)
    body, _ = find_balanced_block(text_nc, open_idx)
    triples = re.findall(r'\(([^)]+)\)', body)
    vals = np.array([[float(x) for x in t.split()] for t in triples])
    assert len(vals) == n
    return vals


def face_area_and_centroid(pts):
    """Newell's method polygon area (magnitude) and centroid, for a planar
    or near-planar polygon given in order."""
    n = len(pts)
    if n < 3:
        return 0.0, pts.mean(axis=0) if n else np.zeros(3)
    c = pts.mean(axis=0)
    normal = np.zeros(3)
    for i in range(n):
        p1 = pts[i]
        p2 = pts[(i + 1) % n]
        normal += np.cross(p1 - c, p2 - c)
    area = 0.5 * np.linalg.norm(normal)
    return area, c


def cell_faces_map(n_cells, faces, owner, neighbour):
    cf = [[] for _ in range(n_cells)]
    n_internal = len(neighbour)
    for fid in range(len(faces)):
        cf[owner[fid]].append(fid)
        if fid < n_internal:
            cf[neighbour[fid]].append(fid)
    return cf


def cell_neighbors(cid, cell_faces, faces, owner, neighbour, n_internal):
    neigh = set()
    for fid in cell_faces[cid]:
        if fid < n_internal:
            other = neighbour[fid] if owner[fid] == cid else owner[fid]
            neigh.add(other)
    return sorted(neigh)


def cell_centroid_and_extent(cid, cell_faces, faces, points):
    pt_ids = set()
    for fid in cell_faces[cid]:
        pt_ids.update(faces[fid])
    pts = points[list(pt_ids)]
    centroid = pts.mean(axis=0)
    x_extent = pts[:, 0].max() - pts[:, 0].min()
    r = np.hypot(pts[:, 1], pts[:, 2])
    r_extent = r.max() - r.min()
    return centroid, x_extent, r_extent, len(pt_ids)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--axis-r-eps", type=float, default=1e-12,
                     help="Points with r <= this are treated as EXACTLY "
                          "on the axis (true r=0).")
    ap.add_argument("--compare-r-tol", type=float, default=0.003,
                     help="Radius tolerance used by the PREVIOUS "
                          "(now-flagged-invalid) extraction, for comparison only.")
    ap.add_argument("--station-x-min", type=float, default=-0.015)
    ap.add_argument("--station-x-max", type=float, default=0.005)
    ap.add_argument("--snapshot", type=str, default=None,
                     help="Optional timestep dir name to pull field data "
                          "(p,rho,T,U) for the proposed axis-cell table.")
    args = ap.parse_args()

    case = Path(args.case_dir)
    poly = case / "constant" / "polyMesh"

    print("Loading mesh geometry...")
    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_labels(poly / "owner")
    neighbour = parse_labels(poly / "neighbour")
    patches = parse_boundary(poly / "boundary")
    n_internal = len(neighbour)
    n_cells = max(owner.max(), neighbour.max()) + 1 if len(neighbour) else owner.max() + 1
    print(f"  n_points={len(points)}  n_faces={len(faces)}  n_cells={n_cells}  "
          f"n_internal_faces={n_internal}\n")

    print("Patches:")
    for name, ptype, nf, sf in patches:
        print(f"  {name:20s} type={ptype:10s} nFaces={nf:8d} startFace={sf}")
    print()

    r_all = np.hypot(points[:, 1], points[:, 2])
    axis_point_mask = r_all <= args.axis_r_eps
    axis_point_ids = set(np.where(axis_point_mask)[0])
    print(f"[1] TRUE AXIS POINTS (r <= {args.axis_r_eps}): {len(axis_point_ids)} of {len(points)} points\n")

    # --- [2] faces containing >=2 axis points (degenerate/axis-adjacent faces) ---
    axis_faces = []
    for fid, f in enumerate(faces):
        n_axis_pts_on_face = sum(1 for p in f if p in axis_point_ids)
        if n_axis_pts_on_face >= 2:
            axis_faces.append((fid, n_axis_pts_on_face))
    print(f"[2] FACES WITH >=2 AXIS POINTS (candidate axis-adjacent faces): {len(axis_faces)}\n")

    axis_cells = set()
    face_area_report = []
    for fid, n_axis_pts in axis_faces:
        pts = points[faces[fid]]
        area, fc = face_area_and_centroid(pts)
        ftype = face_type_for(fid, patches, n_internal)
        axis_cells.add(owner[fid])
        if fid < n_internal:
            axis_cells.add(neighbour[fid])
        face_area_report.append((fid, ftype, area, fc, n_axis_pts))

    print(f"[2b] UNIQUE CELLS ADJACENT TO THESE AXIS FACES: {len(axis_cells)}\n")

    print("Sample of axis-adjacent face details (first 20, sorted by centroid x):")
    face_area_report.sort(key=lambda t: t[3][0])
    print(f"{'faceID':>8} {'type':>20} {'area':>14} {'cx':>10} {'cr':>10} {'n_axis_pts':>11}")
    for fid, ftype, area, fc, nap in face_area_report[:20]:
        cr = np.hypot(fc[1], fc[2])
        print(f"{fid:>8} {ftype:>20} {area:>14.6e} {fc[0]:>10.5f} {cr:>10.6f} {nap:>11}")
    print("  ...")
    for fid, ftype, area, fc, nap in face_area_report[-20:]:
        cr = np.hypot(fc[1], fc[2])
        print(f"{fid:>8} {ftype:>20} {area:>14.6e} {fc[0]:>10.5f} {cr:>10.6f} {nap:>11}")
    print()

    # --- build cell-faces map (only needed for axis cells + their neighbors, to keep it fast) ---
    print("Building cell-face connectivity for axis-adjacent cells and their neighbors...")
    cell_faces = cell_faces_map(n_cells, faces, owner, neighbour)

    axis_cell_info = {}
    for cid in axis_cells:
        centroid, xext, rext, npts = cell_centroid_and_extent(cid, cell_faces, faces, points)
        neigh = cell_neighbors(cid, cell_faces, faces, owner, neighbour, n_internal)
        ftypes = sorted(set(face_type_for(fid, patches, n_internal) for fid in cell_faces[cid]))
        axis_cell_info[cid] = dict(centroid=centroid, x_extent=xext, r_extent=rext,
                                    n_pts=npts, neighbors=neigh, face_types=ftypes,
                                    n_faces=len(cell_faces[cid]))

    print(f"\n[3,4] TRUE AXIS-ADJACENT CELLS: {len(axis_cell_info)}\n")
    print(f"{'cellID':>8} {'cx':>10} {'cr':>10} {'x_ext':>10} {'r_ext':>10} "
          f"{'nFaces':>7} {'nNeigh':>7}  face_types")
    sorted_axis_cells = sorted(axis_cell_info.items(), key=lambda kv: kv[1]["centroid"][0])
    for cid, info in sorted_axis_cells:
        cx, cy, cz = info["centroid"]
        cr = np.hypot(cy, cz)
        print(f"{cid:>8} {cx:>10.5f} {cr:>10.6f} {info['x_extent']:>10.6f} "
              f"{info['r_extent']:>10.6f} {info['n_faces']:>7} {len(info['neighbors']):>7}  "
              f"{info['face_types']}")
    print()

    # --- [3] uniqueness check: is there exactly one axis cell per axial station? ---
    print("[3] UNIQUENESS CHECK: multiple axis-adjacent cells with near-identical x?")
    xs_sorted = sorted(info["centroid"][0] for info in axis_cell_info.values())
    close_pairs = 0
    for i in range(len(xs_sorted) - 1):
        if abs(xs_sorted[i + 1] - xs_sorted[i]) < 1e-6:
            close_pairs += 1
    print(f"  Number of axis-cell x-pairs within 1e-6 m of each other: {close_pairs}")
    print(f"  (0 means the true axis-adjacent cells form a clean monotonic x sequence "
          f"with no duplicates)\n")

    # --- [7,8] station-by-station candidate count comparison ---
    print(f"[7,8] STATION-BY-STATION COMPARISON, x in [{args.station_x_min}, {args.station_x_max}]")
    print("Comparing TRUE axis-adjacent cells vs. previous r<tol selection\n")

    centroids_all = np.zeros((n_cells, 3))
    # compute all cell centroids cheaply via face-average method used previously
    csum = np.zeros((n_cells, 3))
    ccount = np.zeros(n_cells)
    for fid, face in enumerate(faces):
        fc = points[face].mean(axis=0)
        o = owner[fid]
        csum[o] += fc
        ccount[o] += 1
        if fid < n_internal:
            nb = neighbour[fid]
            csum[nb] += fc
            ccount[nb] += 1
    centroids_all = csum / np.maximum(ccount, 1)[:, None]

    r_tol_mask = (np.hypot(centroids_all[:, 1], centroids_all[:, 2]) < args.compare_r_tol) & \
                 (centroids_all[:, 0] > args.station_x_min) & (centroids_all[:, 0] < args.station_x_max)
    r_tol_cells = np.where(r_tol_mask)[0]

    true_axis_in_range = [cid for cid, info in axis_cell_info.items()
                           if args.station_x_min < info["centroid"][0] < args.station_x_max]

    print(f"  TRUE axis-adjacent cells in range: {len(true_axis_in_range)}")
    print(f"  r<{args.compare_r_tol} tolerance-selected cells in range: {len(r_tol_cells)}")
    print(f"  (large disparity indicates the old tolerance-based method was pulling in many "
          f"off-axis cells)\n")

    # bucket both sets into small x-bins to see cells-per-station
    bins = np.linspace(args.station_x_min, args.station_x_max, 21)
    print("  x-bin range          | true_axis_count | r_tol_count")
    for i in range(len(bins) - 1):
        lo, hi = bins[i], bins[i + 1]
        ta = sum(1 for cid in true_axis_in_range if lo <= axis_cell_info[cid]["centroid"][0] < hi)
        rt = sum(1 for cid in r_tol_cells if lo <= centroids_all[cid, 0] < hi)
        print(f"  [{lo:>9.5f}, {hi:>9.5f}) |{ta:>16d} |{rt:>12d}")
    print()

    # --- [9] proposed axis-cell sequence with field data, if snapshot given ---
    if args.snapshot:
        d = case / args.snapshot
        p = parse_scalar_field(d / "p")
        rho = parse_scalar_field(d / "rho")
        T = parse_scalar_field(d / "T")
        U = parse_vector_field(d / "U")
        if p is None or rho is None or T is None or U is None:
            print(f"[9] Could not load fields from snapshot {args.snapshot}")
        else:
            print(f"[9] PROPOSED TRUE-AXIS-CELL SEQUENCE, snapshot t={args.snapshot}")
            gamma, R = 1.4, 287.0
            print(f"{'cellID':>8} {'x':>10} {'r':>10} {'p':>12} {'rho':>12} {'T':>9} "
                  f"{'Ux':>10} {'M':>7}")
            for cid, info in sorted_axis_cells:
                if cid >= len(p):
                    continue
                cx, cy, cz = info["centroid"]
                cr = np.hypot(cy, cz)
                a = np.sqrt(gamma * R * max(T[cid], 1.0))
                M = abs(U[cid, 0]) / a
                print(f"{cid:>8} {cx:>10.5f} {cr:>10.6f} {p[cid]:>12.3f} {rho[cid]:>12.5e} "
                      f"{T[cid]:>9.2f} {U[cid,0]:>10.2f} {M:>7.3f}")
            print()


if __name__ == "__main__":
    main()
