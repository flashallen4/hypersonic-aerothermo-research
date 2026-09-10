#!/usr/bin/env python3
"""
Temporal-independence diagnostic for Project 01 forebody QoIs (read-only,
uses existing saved snapshots only; does not modify mesh/schemes/case or
rerun the solver).

Tracks per snapshot:
  - Stagnation-point wall heat flux: q = kappa(p,T_wall)*(T_wall-T_freestream)/d,
    with kappa via OpenFOAM 11's sutherlandTransport Eucken formula:
    kappa = mu*Cv*(1.32 + 1.77*R/Cv), mu via Sutherland's law.
  - Stagnation-point pressure and wall temperature.
  - Shock stand-off distance: distance from nose tip to the pressure-jump
    location along the stagnation streamline (axis), located as the first
    cell (moving upstream from the nose) where p exceeds a threshold
    fraction of the post-shock stagnation pressure.
  - Forebody surface pressure at several representative wall_cone stations.
  - Base-region breakdown indicators: min p, min rho, count/extent of
    cells with rho < 1% rho_inf (reusing the established threshold).
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


def sutherland_mu(T, As=1.458e-6, Ts=110.4):
    return As * T**1.5 / (T + Ts)


def eucken_kappa(mu, T, R, Cv):
    # kappa = mu*Cv*(1.32 + 1.77*R/Cv)  -- OpenFOAM 11 sutherlandTransportI.H
    return mu * Cv * (1.32 + 1.77 * R / Cv)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--pinf", type=float, default=1197.0)
    ap.add_argument("--rhoinf", type=float, default=0.018407)
    ap.add_argument("--Tinf", type=float, default=226.5)
    ap.add_argument("--Cp", type=float, default=1005.0)
    ap.add_argument("--R", type=float, default=287.0)
    ap.add_argument("--stag-cell", type=int, default=48880)
    ap.add_argument("--stag-face", type=int, default=132600)
    ap.add_argument("--rho-threshold-frac", type=float, default=0.01)
    args = ap.parse_args()

    Cv = args.Cp - args.R

    case = Path(args.case_dir)
    poly = case / "constant" / "polyMesh"

    print("Loading mesh geometry...")
    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_labels(poly / "owner")
    neighbour = parse_labels(poly / "neighbour")
    patches = parse_boundary(poly / "boundary")
    centroids = cell_centroids_corrected(points, faces, owner, neighbour)
    print(f"  n_cells = {len(centroids)}")

    # distance from stagnation cell centroid to stagnation wall face centroid
    stag_cell_cen = centroids[args.stag_cell]
    stag_face_pts = points[faces[args.stag_face]]
    stag_face_cen = stag_face_pts.mean(axis=0)
    d_wall = np.linalg.norm(stag_cell_cen - stag_face_cen)
    print(f"  Stagnation cell {args.stag_cell} centroid: {stag_cell_cen}")
    print(f"  Stagnation wall face {args.stag_face} centroid: {stag_face_cen}")
    print(f"  Wall-normal distance d = {d_wall:.6e} m\n")

    # Forebody pressure stations: pick a few wall_cone faces at different x
    wall_cone = next(p for p in patches if p[0] == "wall_cone")
    _, _, nf_cone, sf_cone = wall_cone
    cone_face_ids = list(range(sf_cone, sf_cone + nf_cone))
    cone_face_x = []
    for fid in cone_face_ids:
        fc = points[faces[fid]].mean(axis=0)
        cone_face_x.append((fid, fc[0]))
    cone_face_x.sort(key=lambda t: t[1])
    # pick faces at ~10%, 30%, 50%, 70%, 90% along the sorted list
    n_stations = 5
    station_faces = [cone_face_x[int(f * (len(cone_face_x) - 1))][0]
                      for f in np.linspace(0.1, 0.9, n_stations)]
    station_x = [dict(cone_face_x).get(fid, None) for fid in station_faces]
    station_cells = [owner[fid] for fid in station_faces]
    print(f"Forebody pressure stations (wall_cone), x-locations: "
          f"{[f'{x:.4f}' for x in station_x]}")
    print(f"  Owner cells: {station_cells}\n")

    # Axis cells for shock stand-off: cells near r=0, x from -0.1 to 0.05 (upstream of nose)
    tdirs = []
    for d in case.iterdir():
        if d.is_dir() and re.match(r'^[\d.eE+-]+$', d.name) and d.name != '0':
            try:
                t = float(d.name)
                tdirs.append((t, d))
            except ValueError:
                continue
    tdirs.sort(key=lambda x: x[0])

    axis_mask = (np.hypot(centroids[:, 1], centroids[:, 2]) < 0.003) & \
                (centroids[:, 0] > -0.20) & (centroids[:, 0] < 0.02)
    axis_idx = np.where(axis_mask)[0]
    axis_x = centroids[axis_idx, 0]
    order = np.argsort(axis_x)
    axis_idx = axis_idx[order]
    axis_x = axis_x[order]
    print(f"Axis cells for shock detection: {len(axis_idx)} cells, "
          f"x range [{axis_x.min():.4f}, {axis_x.max():.4f}]\n")

    rho_thresh = args.rho_threshold_frac * args.rhoinf

    print("=" * 130)
    print("TEMPORAL-INDEPENDENCE TRACKING")
    print("=" * 130)
    header = (f"{'time':>12} {'q_stag[W/m2]':>13} {'p_stag':>9} {'T_wall':>8} | "
              f"{'shock_x':>9} {'standoff':>9} | "
              f"{'p_cone1':>9} {'p_cone2':>9} {'p_cone3':>9} {'p_cone4':>9} {'p_cone5':>9} | "
              f"{'p_min':>9} {'rho_min':>10} {'n_raref':>8}")
    print(header)

    for t, d in tdirs:
        p, p_uni = parse_scalar_field(d / "p")
        T, T_uni = parse_scalar_field(d / "T")
        rho, rho_uni = parse_scalar_field(d / "rho")
        if p is None or T is None or rho is None:
            continue
        n = min(len(p), len(T), len(rho), len(centroids))

        # --- Stagnation point heat flux ---
        T_wall = T[args.stag_cell] if args.stag_cell < n else float('nan')
        p_stag = p[args.stag_cell] if args.stag_cell < n else float('nan')
        mu_wall = sutherland_mu(T_wall)
        kappa_wall = eucken_kappa(mu_wall, T_wall, args.R, Cv)
        dT = T_wall - args.Tinf  # T at wall-adjacent cell minus freestream, as proxy for wall-normal delta
        q_stag = kappa_wall * dT / d_wall

        # --- Forebody pressure stations ---
        p_stations = [p[c] if c < n else float('nan') for c in station_cells]

        # --- Shock stand-off: scan axis cells from far upstream toward nose,
        #     find first cell where p exceeds 2x freestream (shock jump) ---
        shock_x = None
        for idx, xv in zip(axis_idx, axis_x):
            if idx < n and p[idx] > 2.0 * args.pinf:
                shock_x = xv
                break
        standoff = (0.0 - shock_x) if shock_x is not None else float('nan')  # nose at x~0

        # --- Base breakdown indicators ---
        p_min = p[:n].min()
        rho_min = rho[:n].min()
        n_raref = (rho[:n] < rho_thresh).sum()

        line = (f"{t:>12.4e} {q_stag:>13.4e} {p_stag:>9.2f} {T_wall:>8.2f} | "
                f"{str(round(shock_x,5)) if shock_x is not None else 'None':>9} "
                f"{standoff if shock_x is not None else float('nan'):>9.5f} | ")
        line += " ".join(f"{v:>9.2f}" for v in p_stations) + " | "
        line += f"{p_min:>9.4f} {rho_min:>10.4e} {n_raref:>8d}"
        print(line)


if __name__ == "__main__":
    main()
