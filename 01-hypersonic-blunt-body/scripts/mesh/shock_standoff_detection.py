#!/usr/bin/env python3
"""
Robust bow-shock stand-off detection for Project 01 (read-only, uses existing
saved snapshots only; does not modify mesh/schemes/case or rerun the solver).

Method:
  1. PRIMARY: locate the axis cell pair (sorted by x, upstream->downstream)
     with maximum |dp/dx| (one-sided finite difference between consecutive
     axis-cell centroids). This is the standard shock-capturing-scheme
     location criterion for a smeared/captured shock.
  2. CROSS-VALIDATION: repeat using |drho/dx|. Report agreement/disagreement
     with the pressure-gradient location (in both x-distance and cell-index
     terms).
  3. RANKINE-HUGONIOT SANITY CHECK: using a representative freestream
     (pre-shock) axis cell far upstream and a representative post-shock axis
     cell just downstream of the detected shock location, compute observed
     p2/p1, rho2/rho1, T2/T1 and compare against the theoretical normal-shock
     relations for the specified freestream Mach number and gamma.

Explicitly REJECTED prior method (documented, not used here): first axis
cell where p > 2*p_inf scanning from upstream. That method is NOT invoked
by this script.
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


def normal_shock_ratios(M1, gamma):
    p2_p1 = (2 * gamma * M1**2 - (gamma - 1)) / (gamma + 1)
    rho2_rho1 = ((gamma + 1) * M1**2) / ((gamma - 1) * M1**2 + 2)
    T2_T1 = p2_p1 / rho2_rho1
    return p2_p1, rho2_rho1, T2_T1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--pinf", type=float, default=1197.0)
    ap.add_argument("--rhoinf", type=float, default=0.018407)
    ap.add_argument("--Tinf", type=float, default=226.5)
    ap.add_argument("--Minf", type=float, default=7.0)
    ap.add_argument("--gamma", type=float, default=1.4)
    ap.add_argument("--R", type=float, default=287.0)
    ap.add_argument("--axis-r-tol", type=float, default=0.003)
    ap.add_argument("--axis-x-min", type=float, default=-0.20)
    ap.add_argument("--axis-x-max", type=float, default=0.02)
    ap.add_argument("--freestream-x-max", type=float, default=-0.12,
                     help="Axis cells with x < this are treated as the "
                          "representative pre-shock/freestream sample.")
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

    axis_mask = (np.hypot(centroids[:, 1], centroids[:, 2]) < args.axis_r_tol) & \
                (centroids[:, 0] > args.axis_x_min) & (centroids[:, 0] < args.axis_x_max)
    axis_idx = np.where(axis_mask)[0]
    axis_x = centroids[axis_idx, 0]
    order = np.argsort(axis_x)
    axis_idx = axis_idx[order]
    axis_x = axis_x[order]
    print(f"Axis cells: {len(axis_idx)}, x range [{axis_x.min():.5f}, {axis_x.max():.5f}]\n")

    p2_p1_theory, rho2_rho1_theory, T2_T1_theory = normal_shock_ratios(args.Minf, args.gamma)
    print(f"Theoretical normal-shock ratios at M_inf={args.Minf}, gamma={args.gamma}:")
    print(f"  p2/p1   = {p2_p1_theory:.4f}")
    print(f"  rho2/rho1 = {rho2_rho1_theory:.4f}")
    print(f"  T2/T1   = {T2_T1_theory:.4f}\n")

    freestream_sel = axis_x < args.freestream_x_max
    if not freestream_sel.any():
        print("WARNING: no axis cells satisfy freestream-x-max criterion; "
              "widen --axis-x-min or relax --freestream-x-max.")
        return

    tdirs = []
    for d in case.iterdir():
        if d.is_dir() and re.match(r'^[\d.eE+-]+$', d.name) and d.name != '0':
            try:
                t = float(d.name)
                tdirs.append((t, d))
            except ValueError:
                continue
    tdirs.sort(key=lambda x: x[0])

    print("=" * 150)
    print("SHOCK STAND-OFF DETECTION (peak |dp/dx| primary, peak |drho/dx| cross-check)")
    print("=" * 150)
    header = (f"{'time':>12} | {'x_shock_p':>10} {'standoff_p':>10} | "
              f"{'x_shock_rho':>11} {'standoff_rho':>12} {'agree(cells)':>12} | "
              f"{'p2/p1_obs':>9} {'rho2/rho1_obs':>13} {'T2/T1_obs':>9}")
    print(header)

    for t, d in tdirs:
        p = parse_scalar_field(d / "p")
        rho = parse_scalar_field(d / "rho")
        T = parse_scalar_field(d / "T")
        if p is None or rho is None or T is None:
            continue
        n = min(len(p), len(rho), len(T), len(centroids))
        valid = axis_idx < n
        aidx = axis_idx[valid]
        ax = axis_x[valid]
        if len(aidx) < 3:
            continue

        p_axis = p[aidx]
        rho_axis = rho[aidx]
        T_axis = T[aidx]

        dx = np.diff(ax)
        dp = np.diff(p_axis)
        drho = np.diff(rho_axis)
        dpdx = dp / dx
        drhodx = drho / dx

        i_p = np.argmax(np.abs(dpdx))
        i_rho = np.argmax(np.abs(drhodx))

        x_shock_p = 0.5 * (ax[i_p] + ax[i_p + 1])
        x_shock_rho = 0.5 * (ax[i_rho] + ax[i_rho + 1])
        standoff_p = -x_shock_p
        standoff_rho = -x_shock_rho
        cell_agree = i_rho - i_p

        fs_local = ax < args.freestream_x_max
        if not fs_local[:len(aidx)].any():
            p1 = args.pinf
            rho1 = args.rhoinf
        else:
            p1 = p_axis[fs_local[:len(aidx)]].mean()
            rho1 = rho_axis[fs_local[:len(aidx)]].mean()

        post_idx = min(i_p + 2, len(aidx) - 1)
        p2 = p_axis[post_idx]
        rho2 = rho_axis[post_idx]
        T1 = args.Tinf
        T2 = T_axis[post_idx]

        p2_p1_obs = p2 / p1
        rho2_rho1_obs = rho2 / rho1
        T2_T1_obs = T2 / T1

        line = (f"{t:>12.4e} | {x_shock_p:>10.5f} {standoff_p:>10.5f} | "
                f"{x_shock_rho:>11.5f} {standoff_rho:>12.5f} {cell_agree:>12d} | "
                f"{p2_p1_obs:>9.3f} {rho2_rho1_obs:>13.3f} {T2_T1_obs:>9.3f}")
        print(line)


if __name__ == "__main__":
    main()
