#!/usr/bin/env python3
"""
Read-only local neighborhood diagnostic around specific cell IDs.
Reports cell geometry (centroid, volume, face areas, extent), field values
(T, p, rho, U) at a given timestep, neighbor connectivity, and wall-patch
adjacency, for a BFS neighborhood around seed cells.
Does not modify mesh, case, or run the solver.
"""

import re
import argparse
from pathlib import Path
import numpy as np


def find_balanced_block(text, start_idx):
    assert text[start_idx] == '('
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


def parse_scalar_field(path):
    text = Path(path).read_text()
    text_nc = re.sub(r'//.*', '', text)
    m = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\n?\s*(\d+)\s*\n', text_nc)
    if m is None:
        m2 = re.search(r'internalField\s+uniform\s+([\-\d.eE]+)', text_nc)
        return (None, float(m2.group(1))) if m2 else (None, None)
    n = int(m.group(1))
    open_idx = text_nc.index('(', m.end() - 1)
    body, _ = find_balanced_block(text_nc, open_idx)
    vals = np.array([float(x) for x in body.split()])
    assert len(vals) == n
    return vals, None


def parse_vector_field(path):
    text = Path(path).read_text()
    text_nc = re.sub(r'//.*', '', text)
    m = re.search(r'internalField\s+nonuniform\s+List<vector>\s*\n?\s*(\d+)\s*\n', text_nc)
    if m is None:
        return None
    n = int(m.group(1))
    open_idx = text_nc.index('(', m.end() - 1)
    body, _ = find_balanced_block(text_nc, open_idx)
    tuples = re.findall(r'\(([^()]+)\)', body)
    vecs = np.array([[float(x) for x in t.split()] for t in tuples])
    assert len(vecs) == n
    return vecs


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


def face_area_normal(pts):
    """Newell's method: returns (area, unit_normal, centroid)."""
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("timestep")
    ap.add_argument("--seed-cells", type=int, nargs="+", default=[4632, 4633])
    ap.add_argument("--bfs-depth", type=int, default=2)
    args = ap.parse_args()

    case = Path(args.case_dir)
    poly = case / "constant" / "polyMesh"
    tdir = case / args.timestep

    print(f"Case: {case}, timestep: {args.timestep}")
    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_labels(poly / "owner")
    neighbour = parse_labels(poly / "neighbour")
    patches = parse_boundary(poly / "boundary")
    n_internal = len(neighbour)

    T, T_uni = parse_scalar_field(tdir / "T")
    p, p_uni = parse_scalar_field(tdir / "p")
    rho, rho_uni = parse_scalar_field(tdir / "rho")
    U = parse_vector_field(tdir / "U")

    def field_val(arr, uni, cid):
        if arr is not None:
            return arr[cid]
        return uni

    # Build cell -> faces map (once)
    print("Building cell->face map (may take a moment)...")
    cell_faces = {}
    for fid, o in enumerate(owner):
        cell_faces.setdefault(o, []).append(fid)
    for fid, n in enumerate(neighbour):
        cell_faces.setdefault(n, []).append(fid)

    def cell_neighbors(cid):
        result = []
        for fid in cell_faces.get(cid, []):
            if fid < n_internal:
                other = neighbour[fid] if owner[fid] == cid else owner[fid]
                result.append((other, fid))
        return result

    def cell_geom(cid):
        fids = cell_faces.get(cid, [])
        vert_ids = set()
        for fid in fids:
            vert_ids.update(faces[fid])
        pts = points[list(vert_ids)]
        centroid = pts.mean(axis=0)
        bmin, bmax = pts.min(axis=0), pts.max(axis=0)
        # volume via divergence theorem
        vol = 0.0
        for fid in fids:
            fpts = points[faces[fid]]
            area, unit_n, fcen = face_area_normal(fpts)
            sign = 1.0 if owner[fid] == cid else -1.0
            vol += sign * np.dot(fcen, unit_n) * area
        vol = vol / 3.0
        return centroid, bmin, bmax, vol, fids

    # BFS neighborhood
    visited = set(args.seed_cells)
    frontier = set(args.seed_cells)
    for depth in range(args.bfs_depth):
        new_frontier = set()
        for cid in frontier:
            for other, fid in cell_neighbors(cid):
                if other not in visited:
                    new_frontier.add(other)
        visited.update(new_frontier)
        frontier = new_frontier

    print(f"\nNeighborhood size (BFS depth {args.bfs_depth} from seeds {args.seed_cells}): {len(visited)} cells")

    # Report seeds in detail first
    for cid in args.seed_cells:
        centroid, bmin, bmax, vol, fids = cell_geom(cid)
        x, y, z = centroid
        r = np.hypot(y, z)
        print(f"\n{'='*70}")
        print(f"SEED CELL {cid}: centroid x={x:.6f} r={r:.6f}")
        print(f"  bbox extent: {bmax - bmin}")
        print(f"  volume: {vol:.6e} m^3")
        print(f"  fields: T={field_val(T,T_uni,cid):.4f} K  p={field_val(p,p_uni,cid):.4f} Pa  "
              f"rho={field_val(rho,rho_uni,cid):.6e}  U=({U[cid,0]:.3f},{U[cid,1]:.3f},{U[cid,2]:.3f})")
        print(f"  faces ({len(fids)}):")
        for fid in fids:
            fpts = points[faces[fid]]
            area, unit_n, fcen = face_area_normal(fpts)
            cls, pname = classify_face(fid, n_internal, patches)
            if cls == "internal":
                other = neighbour[fid] if owner[fid] == cid else owner[fid]
                print(f"    face {fid}: area={area:.4e}  {cls:10s} -> neighbor cell {other} "
                      f"(T={field_val(T,T_uni,other):.2f} p={field_val(p,p_uni,other):.2f} "
                      f"rho={field_val(rho,rho_uni,other):.4e} "
                      f"U=({U[other,0]:.2f},{U[other,1]:.2f},{U[other,2]:.2f}))")
            else:
                print(f"    face {fid}: area={area:.4e}  boundary:{pname}")

    # Report full neighborhood table sorted by x then r
    print(f"\n{'='*70}")
    print(f"FULL NEIGHBORHOOD TABLE ({len(visited)} cells)")
    print(f"{'='*70}")
    rows = []
    for cid in visited:
        centroid, bmin, bmax, vol, fids = cell_geom(cid)
        x, y, z = centroid
        r = np.hypot(y, z)
        ext = bmax - bmin
        nz = [e for e in ext if e > 1e-12]
        ar = max(nz) / min(nz) if nz else float('nan')
        # wall adjacency
        wall_patches = set()
        for fid in fids:
            cls, pname = classify_face(fid, n_internal, patches)
            if cls == "boundary" and pname in ("wall_cone", "wall_base", "wall_nose"):
                wall_patches.add(pname)
        rows.append((cid, x, r, vol, ar,
                     field_val(T, T_uni, cid), field_val(p, p_uni, cid),
                     field_val(rho, rho_uni, cid), U[cid, 0],
                     np.linalg.norm(U[cid]), ",".join(wall_patches) if wall_patches else "-"))

    rows.sort(key=lambda t: (t[1], t[2]))
    print(f"{'cell':>6} {'x':>10} {'r':>10} {'vol':>12} {'AR':>8} "
          f"{'T':>9} {'p':>10} {'rho':>12} {'Ux':>10} {'|U|':>10} {'wall'}")
    for cid, x, r, vol, ar, Tv, pv, rhov, ux, umag, wall in rows:
        print(f"{cid:>6} {x:>10.6f} {r:>10.6f} {vol:>12.4e} {ar:>8.1f} "
              f"{Tv:>9.2f} {pv:>10.2f} {rhov:>12.4e} {ux:>10.2f} {umag:>10.2f} {wall}")

    # Reference comparison: a "normal" BL cell on the cone wall away from the corner
    print(f"\n{'='*70}")
    print("REFERENCE: normal cone-wall BL cell far from corner (for contrast)")
    print(f"{'='*70}")
    # find a cell whose centroid x is mid-cone (~0.2) and small r extent (BL-like)
    target_x = 0.2
    best_cid, best_dist = None, 1e9
    for cid in range(len(owner) if False else 0):
        pass
    # scan a sample of cells via faces near mid-cone wall patch
    wall_cone_patch = next((pp for pp in patches if pp[0] == "wall_cone"), None)
    if wall_cone_patch:
        name, ptype, nfaces, startface = wall_cone_patch
        # pick a face roughly in the middle of the patch
        mid_fid = startface + nfaces // 2
        fpts = points[faces[mid_fid]]
        _, _, fcen = face_area_normal(fpts)
        ref_cid = owner[mid_fid]
        centroid, bmin, bmax, vol, fids = cell_geom(ref_cid)
        x, y, z = centroid
        r = np.hypot(y, z)
        ext = bmax - bmin
        nz = [e for e in ext if e > 1e-12]
        ar = max(nz) / min(nz) if nz else float('nan')
        print(f"Reference cell {ref_cid}: x={x:.6f} r={r:.6f} vol={vol:.4e} AR={ar:.1f} "
              f"T={field_val(T,T_uni,ref_cid):.2f} p={field_val(p,p_uni,ref_cid):.2f} "
              f"rho={field_val(rho,rho_uni,ref_cid):.4e} "
              f"U=({U[ref_cid,0]:.2f},{U[ref_cid,1]:.2f},{U[ref_cid,2]:.2f})")
    else:
        print("wall_cone patch not found.")


if __name__ == "__main__":
    main()
