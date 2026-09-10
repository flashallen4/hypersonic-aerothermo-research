#!/usr/bin/env python3
"""
Read-only: estimate the gradient-length-local Knudsen number (Boyd et al.
criterion) in and around the rarefied region, to assess continuum validity.
Does not modify mesh, schemes, or run the solver.

Kn_GLL = (lambda / Q) * |dQ/dx|_max_over_neighbors

lambda (mean free path) computed via a standard kinetic-theory relation
using local rho, T, and Sutherland-law viscosity (mu).

Breakdown criteria (Boyd et al. 1995, widely used):
  Kn_GLL > 0.05  -> continuum breakdown likely in that region
  Kn_GLL > 0.1   -> Navier-Stokes clearly unreliable
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


def build_neighbor_map(owner, neighbour):
    n_internal = len(neighbour)
    nbmap = {}
    for fid in range(n_internal):
        o, nb = owner[fid], neighbour[fid]
        nbmap.setdefault(o, []).append(nb)
        nbmap.setdefault(nb, []).append(o)
    return nbmap


def sutherland_mu(T, As=1.458e-6, Ts=110.4):
    return As * T**1.5 / (T + Ts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--timesteps", nargs="+", required=True,
                     help="Specific timestep directory names to analyze")
    ap.add_argument("--R", type=float, default=287.0)
    ap.add_argument("--rho-threshold-frac", type=float, default=0.01)
    ap.add_argument("--rhoinf", type=float, default=0.018407)
    args = ap.parse_args()

    case = Path(args.case_dir)
    poly = case / "constant" / "polyMesh"

    print("Loading mesh geometry...")
    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_labels(poly / "owner")
    neighbour = parse_labels(poly / "neighbour")
    centroids = cell_centroids_corrected(points, faces, owner, neighbour)
    nbmap = build_neighbor_map(owner, neighbour)
    print(f"  n_cells = {len(centroids)}\n")

    rho_thresh = args.rho_threshold_frac * args.rhoinf

    for tname in args.timesteps:
        d = case / tname
        rho, _ = parse_scalar_field(d / "rho")
        T, _ = parse_scalar_field(d / "T")
        if rho is None or T is None:
            print(f"t={tname}: missing fields, skipping")
            continue

        n = min(len(rho), len(T), len(centroids))
        rho, T, cen = rho[:n], T[:n], centroids[:n]

        mu = sutherland_mu(T)
        # mean free path via standard kinetic-theory relation:
        # lambda = mu / rho * sqrt(pi / (2 * R * T))
        lam = (mu / rho) * np.sqrt(np.pi / (2 * args.R * T))

        # gradient-length-local Knudsen number based on density gradient
        Kn = np.zeros(n)
        for cid in range(n):
            if cid not in nbmap:
                continue
            neighbors = nbmap[cid]
            max_grad_term = 0.0
            for nb in neighbors:
                if nb >= n:
                    continue
                dQ = abs(rho[cid] - rho[nb])
                dx = np.linalg.norm(cen[cid] - cen[nb])
                if dx > 0 and rho[cid] > 0:
                    grad_term = (lam[cid] / rho[cid]) * (dQ / dx)
                    max_grad_term = max(max_grad_term, grad_term)
            Kn[cid] = max_grad_term

        mask = rho < rho_thresh
        n_rarefied = mask.sum()

        print(f"=== t={tname} ===")
        print(f"  Rarefied cells (rho<{rho_thresh:.3e}): {n_rarefied}")
        if n_rarefied > 0:
            Kn_rarefied = Kn[mask]
            print(f"  Kn_GLL in rarefied region: min={Kn_rarefied.min():.4f} "
                  f"mean={Kn_rarefied.mean():.4f} max={Kn_rarefied.max():.4f}")
            print(f"  Mean free path (lambda) in rarefied region: "
                  f"min={lam[mask].min():.4e} m, max={lam[mask].max():.4e} m")
            frac_breakdown = (Kn_rarefied > 0.05).sum() / len(Kn_rarefied) * 100
            frac_severe = (Kn_rarefied > 0.1).sum() / len(Kn_rarefied) * 100
            print(f"  Fraction of rarefied cells with Kn>0.05 (breakdown likely): {frac_breakdown:.1f}%")
            print(f"  Fraction of rarefied cells with Kn>0.1 (NS clearly unreliable): {frac_severe:.1f}%")
        print(f"  Global Kn_GLL: min={Kn.min():.4f} mean={Kn.mean():.4f} max={Kn.max():.4f}")
        print()


if __name__ == "__main__":
    main()
