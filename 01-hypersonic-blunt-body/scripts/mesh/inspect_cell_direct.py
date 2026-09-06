#!/usr/bin/env python3
"""
Diagnostic (read-only): directly inspect a specific cell by ID in one case,
and separately find the geometrically closest cell to a given (x,r) target
in another case for comparison. No mesh modification.

Usage:
    python3 inspect_cell_direct.py <case_dir> --cell-id 4633
    python3 inspect_cell_direct.py <case_dir> --x 0.416627 --r 0.149440 --tol 0.00005
"""
import re, argparse
from pathlib import Path
import numpy as np


def extract_block(text):
    idx = text.index('}')
    text = text[idx+1:]
    start = text.index('(')
    depth = 0
    for i in range(start, len(text)):
        if text[i] == '(':
            depth += 1
        elif text[i] == ')':
            depth -= 1
            if depth == 0:
                return text[start+1:i]
    raise ValueError


def read_points(p):
    b = extract_block(Path(p).read_text())
    return np.array([[float(v) for v in e.split()] for e in re.findall(r'\(([^()]+)\)', b)])


def read_labels(p):
    b = extract_block(Path(p).read_text())
    return np.array([int(n) for n in re.findall(r'-?\d+', b)])


def read_faces(p):
    b = extract_block(Path(p).read_text())
    return [[int(v) for v in m.group(2).split()] for m in re.finditer(r'(\d+)\(([^()]+)\)', b)]


def read_boundary(p):
    b = extract_block(Path(p).read_text())
    out = []
    for m in re.finditer(r'(\w+)\s*\{[^{}]*?type\s+(\w+);[^{}]*?nFaces\s+(\d+);[^{}]*?startFace\s+(\d+);[^{}]*?\}', b, re.S):
        n, t, nf, sf = m.groups()
        out.append((n, t, int(nf), int(sf)))
    return out


def classify(fid, n_int, patches):
    if fid < n_int:
        return "internal"
    for n, t, nf, sf in patches:
        if sf <= fid < sf + nf:
            return f"boundary:{n}"
    return "boundary:UNKNOWN"


def area(pts):
    c = pts.mean(axis=0)
    nv = np.zeros(3)
    for i in range(len(pts)):
        nv += np.cross(pts[i]-c, pts[(i+1) % len(pts)]-c)
    return 0.5*np.linalg.norm(nv)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--cell-id", type=int, default=None)
    ap.add_argument("--x", type=float, default=None)
    ap.add_argument("--r", type=float, default=None)
    ap.add_argument("--tol", type=float, default=0.00005)
    args = ap.parse_args()

    pm = Path(args.case_dir) / "constant" / "polyMesh"
    points = read_points(pm / "points")
    faces = read_faces(pm / "faces")
    owner = read_labels(pm / "owner")
    neighbour = read_labels(pm / "neighbour")
    patches = read_boundary(pm / "boundary")
    n_int = len(neighbour)

    def faces_of(cid):
        return [f for f, o in enumerate(owner) if o == cid] + \
               [f for f, n in enumerate(neighbour) if n == cid]

    def report(cid):
        fids = faces_of(cid)
        vids = set()
        for f in fids:
            vids.update(faces[f])
        pts = points[list(vids)]
        cen = pts.mean(axis=0)
        bmin, bmax = pts.min(axis=0), pts.max(axis=0)
        x, y, z = cen
        r = np.hypot(y, z)
        print(f"\nCell {cid}: centroid x={x:.6f} r={r:.6f}")
        print(f"  bbox min={bmin}  max={bmax}  extent={bmax-bmin}")
        print(f"  faces: {len(fids)}")
        for f in fids:
            v = faces[f]
            p = points[v]
            a = area(p)
            cls = classify(f, n_int, patches)
            fx = p[:, 0].mean(); fr = np.hypot(p[:, 1], p[:, 2]).mean()
            nb = None
            if f < n_int:
                nb = neighbour[f] if owner[f] == cid else owner[f]
            print(f"    face {f}: area={a:.6e} x={fx:.6f} r={fr:.6f} {cls:25s} neighbor={nb}")
        return cen

    if args.cell_id is not None:
        report(args.cell_id)
    elif args.x is not None:
        # find cell(s) whose bbox contains (x,r) within tol
        candidates = []
        n_cells = max(owner.max(), neighbour.max()) + 1
        # crude but exhaustive-safe: check via faces near target first
        near_faces = []
        for fid, verts in enumerate(faces):
            p = points[verts]
            r_ = np.hypot(p[:, 1], p[:, 2])
            x_ = p[:, 0]
            if np.any((np.abs(x_-args.x) < args.tol) & (np.abs(r_-args.r) < args.tol)):
                near_faces.append(fid)
        cells = set()
        for fid in near_faces:
            cells.add(owner[fid])
            if fid < n_int:
                cells.add(neighbour[fid])
        print(f"Found {len(cells)} candidate cells within tol={args.tol}")
        results = []
        for cid in cells:
            fids = faces_of(cid)
            vids = set()
            for f in fids:
                vids.update(faces[f])
            pts = points[list(vids)]
            cen = pts.mean(axis=0)
            ext = pts.max(axis=0) - pts.min(axis=0)
            results.append((cid, cen, ext))
        results.sort(key=lambda t: min([e for e in t[2] if e > 1e-12]))
        print("Sorted by smallest nonzero extent (thinnest first):")
        for cid, cen, ext in results[:10]:
            nz = [e for e in ext if e > 1e-12]
            ar = max(nz)/min(nz) if nz else float('nan')
            print(f"  Cell {cid}: x={cen[0]:.6f} r={np.hypot(cen[1],cen[2]):.6f} extent={ext} AR={ar:.2f}")
        if results:
            report(results[0][0])


if __name__ == "__main__":
    main()
