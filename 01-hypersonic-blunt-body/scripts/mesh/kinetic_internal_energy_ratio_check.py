#!/usr/bin/env python3
"""
Read-only: compute the ratio of kinetic energy to sensible internal energy
per unit mass (KE/IE = 0.5*|U|^2 / (Cv*T)) at the pressure-minimum cell and
its neighborhood, across snapshots. A large and growing ratio indicates
that deriving T from conserved total energy involves subtracting two large,
nearly-equal quantities (total energy, kinetic energy) to recover a small
remainder (sensible internal energy) -- a classic catastrophic-cancellation
failure mode in density-based compressible solvers as rho->0.

Does not modify mesh, schemes, fields, or run the solver.
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
    ap.add_argument("--Cv", type=float, default=718.0)  # Cp - R = 1005 - 287
    args = ap.parse_args()

    case = Path(args.case_dir)
    poly = case / "constant" / "polyMesh"

    print("Loading mesh geometry...")
    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_labels(poly / "owner")
    neighbour = parse_labels(poly / "neighbour")
    centroids = cell_centroids_corrected(points, faces, owner, neighbour)
    print(f"  n_cells = {len(centroids)}\n")

    tdirs = []
    for d in case.iterdir():
        if d.is_dir() and re.match(r'^[\d.eE+-]+$', d.name) and d.name != '0':
            try:
                t = float(d.name)
                tdirs.append((t, d))
            except ValueError:
                continue
    tdirs.sort(key=lambda x: x[0])

    print(f"{'time':>12} {'p_min_cell':>10} {'T':>10} {'|U|':>10} {'KE/mass':>12} "
          f"{'IE/mass(CvT)':>12} {'KE/IE ratio':>12} {'p':>10}")
    print("-" * 100)

    for t, d in tdirs:
        p, p_uni = parse_scalar_field(d / "p")
        T, T_uni = parse_scalar_field(d / "T")
        U = parse_vector_field(d / "U")
        if p is None or T is None or U is None:
            continue
        n = min(len(p), len(T), len(U))
        p, T, U = p[:n], T[:n], U[:n]

        idx_pmin = np.argmin(p)
        Umag = np.linalg.norm(U[idx_pmin])
        KE = 0.5 * Umag**2
        IE = args.Cv * T[idx_pmin]
        ratio = KE / IE if IE > 0 else float('inf')

        print(f"{t:>12.4e} {idx_pmin:>10d} {T[idx_pmin]:>10.3f} {Umag:>10.3f} "
              f"{KE:>12.4e} {IE:>12.4e} {ratio:>12.4f} {p[idx_pmin]:>10.4f}")

    print("\n--- Also checking global max KE/IE ratio (may not be at p_min cell) ---")
    for t, d in tdirs:
        p, p_uni = parse_scalar_field(d / "p")
        T, T_uni = parse_scalar_field(d / "T")
        U = parse_vector_field(d / "U")
        if p is None or T is None or U is None:
            continue
        n = min(len(p), len(T), len(U))
        T, U = T[:n], U[:n]
        Umag = np.linalg.norm(U, axis=1)
        KE = 0.5 * Umag**2
        IE = args.Cv * np.maximum(T, 1e-6)
        ratio = KE / IE
        idx_max = np.argmax(ratio)
        print(f"t={t:.4e}: max KE/IE ratio = {ratio[idx_max]:.4f} at cell {idx_max} "
              f"(T={T[idx_max]:.2f}, |U|={Umag[idx_max]:.2f})")


if __name__ == "__main__":
    main()
