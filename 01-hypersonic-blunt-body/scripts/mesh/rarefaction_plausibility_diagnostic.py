#!/usr/bin/env python3
"""
Read-only diagnostic: physical plausibility and spatial/temporal progression
of the low-pressure/low-density region near the cone/base corner and wake,
across all saved snapshots of a case. Does not modify mesh, schemes, or
run the solver.

Part 1: global min/max of p, p/p_inf, rho, rho/rho_inf, T, T/T_inf, |U|, Ux,
        local Mach number; location of minima; ideal-gas consistency check
        (p vs rho*R*T residual).
Part 2: for each snapshot, spatial extent (x,r bounding box, cell count) of
        cells below a low-density threshold, tracked over time, to assess
        whether the region is bounded/spreading/strengthening, and whether
        it is approaching the outlet.
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


def parse_scalar_field(path):
    if not path.exists():
        return None, None
    text = path.read_text()
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


def cell_centroids_corrected(points, faces, owner, neighbour):
    n_cells = max(owner.max(), neighbour.max()) + 1 if len(neighbour) else owner.max() + 1
    csum = np.zeros((n_cells, 3))
    ccount = np.zeros(n_cells)
    n_internal = len(neighbour)
    for fid, face in enumerate(faces):
        fc = points[face].mean(axis=0)
        o = owner[fid]
        csum[o] += fc
        ccount[o] += 1
        if fid < n_internal:
            nb = neighbour[fid]
            csum[nb] += fc
            ccount[nb] += 1
    return csum / np.maximum(ccount, 1)[:, None]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--pinf", type=float, default=1197.0)
    ap.add_argument("--rhoinf", type=float, default=0.018407)
    ap.add_argument("--Tinf", type=float, default=226.5)
    ap.add_argument("--Uinf", type=float, default=2111.72)
    ap.add_argument("--gamma", type=float, default=1.4)
    ap.add_argument("--R", type=float, default=287.0)
    ap.add_argument("--rho-threshold-frac", type=float, default=0.01,
                     help="fraction of rho_inf below which a cell counts as 'rarefied'")
    args = ap.parse_args()

    case = Path(args.case_dir)
    poly = case / "constant" / "polyMesh"

    print("Loading mesh geometry...")
    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_labels(poly / "owner")
    neighbour = parse_labels(poly / "neighbour")
    centroids = cell_centroids_corrected(points, faces, owner, neighbour)
    print(f"  n_cells = {len(centroids)}")

    tdirs = []
    for d in case.iterdir():
        if d.is_dir() and re.match(r'^[\d.eE+-]+$', d.name) and d.name != '0':
            try:
                t = float(d.name)
                tdirs.append((t, d))
            except ValueError:
                continue
    tdirs.sort(key=lambda x: x[0])
    print(f"  {len(tdirs)} non-zero timestep directories found\n")

    R = args.R
    rho_thresh = args.rho_threshold_frac * args.rhoinf

    print("=" * 100)
    print("PART 1: PHYSICAL PLAUSIBILITY (global extrema per snapshot)")
    print("=" * 100)
    header = (f"{'time':>12} {'p_min':>10} {'p/pinf_min':>11} {'rho_min':>12} "
              f"{'rho/rhoinf':>11} {'T_at_pmin':>10} {'|U|_max':>9} {'M_max':>7} "
              f"{'loc(x,r)':>18} {'gaslaw_err%':>11}")
    print(header)

    part1_rows = []
    for t, d in tdirs:
        p, p_uni = parse_scalar_field(d / "p")
        T, T_uni = parse_scalar_field(d / "T")
        rho, rho_uni = parse_scalar_field(d / "rho")
        U = parse_vector_field(d / "U")
        if p is None or T is None or rho is None or U is None:
            continue
        n = min(len(p), len(T), len(rho), len(U), len(centroids))
        p, T, rho, U = p[:n], T[:n], rho[:n], U[:n]
        cen = centroids[:n]

        Umag = np.linalg.norm(U, axis=1)
        a = np.sqrt(args.gamma * R * np.maximum(T, 1e-6))
        Mach = Umag / a

        idx_pmin = np.argmin(p)
        x_pmin, y_pmin, z_pmin = cen[idx_pmin]
        r_pmin = np.hypot(y_pmin, z_pmin)

        # ideal gas law consistency: p_calc = rho*R*T vs actual p
        p_calc = rho * R * T
        gaslaw_err = np.abs(p - p_calc) / np.maximum(np.abs(p), 1e-30)
        gaslaw_err_at_pmin = gaslaw_err[idx_pmin] * 100

        row = dict(t=t, p_min=p.min(), rho_min=rho.min(), T_at_pmin=T[idx_pmin],
                   Umag_max=Umag.max(), Mach_max=Mach.max(),
                   x_pmin=x_pmin, r_pmin=r_pmin, gaslaw_err=gaslaw_err_at_pmin,
                   p=p, rho=rho, T=T, Umag=Umag, Mach=Mach, cen=cen)
        part1_rows.append(row)

        print(f"{t:>12.4e} {p.min():>10.4f} {p.min()/args.pinf:>11.4e} {rho.min():>12.4e} "
              f"{rho.min()/args.rhoinf:>11.4e} {T[idx_pmin]:>10.3f} {Umag.max():>9.2f} "
              f"{Mach.max():>7.2f} ({x_pmin:.5f},{r_pmin:.5f}) {gaslaw_err_at_pmin:>10.4f}")

    print("\n--- Rate of decrease: p_min and rho_min ratio between successive snapshots ---")
    for i in range(1, len(part1_rows)):
        prev, cur = part1_rows[i - 1], part1_rows[i]
        dt = cur['t'] - prev['t']
        if prev['p_min'] > 0 and dt > 0:
            p_ratio = cur['p_min'] / prev['p_min']
            rho_ratio = cur['rho_min'] / prev['rho_min']
            print(f"  t={prev['t']:.3e}->{cur['t']:.3e} (dt={dt:.3e}): "
                  f"p_min ratio={p_ratio:.4f}  rho_min ratio={rho_ratio:.4f}")

    print("\n" + "=" * 100)
    print(f"PART 2: SPATIAL PROGRESSION of rarefied region (rho < {args.rho_threshold_frac}*rho_inf = {rho_thresh:.4e})")
    print("=" * 100)
    header2 = (f"{'time':>12} {'n_cells':>8} {'x_min':>10} {'x_max':>10} "
               f"{'r_min':>10} {'r_max':>10} {'x_extent':>10} {'r_extent':>10}")
    print(header2)
    for row in part1_rows:
        mask = row['rho'] < rho_thresh
        n_rarefied = mask.sum()
        if n_rarefied == 0:
            print(f"{row['t']:>12.4e} {0:>8} {'--':>10} {'--':>10} {'--':>10} {'--':>10} {'--':>10} {'--':>10}")
            continue
        cen_r = row['cen'][mask]
        x_r = cen_r[:, 0]
        r_r = np.hypot(cen_r[:, 1], cen_r[:, 2])
        print(f"{row['t']:>12.4e} {n_rarefied:>8} {x_r.min():>10.5f} {x_r.max():>10.5f} "
              f"{r_r.min():>10.5f} {r_r.max():>10.5f} {x_r.max()-x_r.min():>10.5f} {r_r.max()-r_r.min():>10.5f}")

    print("\nDone.")


if __name__ == "__main__":
    main()
