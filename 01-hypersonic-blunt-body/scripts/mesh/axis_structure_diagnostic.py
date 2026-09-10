#!/usr/bin/env python3
"""
Axis flow-structure diagnostic for Project 01 (read-only, existing saved
snapshots only; no case modification, no solver rerun).

Purpose: understand the actual axial solution structure BEFORE defining any
automated shock-location criterion. Does NOT declare a shock location. Does
NOT apply arbitrary absolute thresholds (no 2x p_inf, no 10%-80% rise, no
fixed wall-distance cutoff, no uniform-grid interpolation).

For each snapshot, computes the full stagnation-axis profile (x, p, rho, T,
M, |dp/dx|, |drho/dx|) and finds ALL local maxima (not just the global max)
in the two gradient series. Each local-maximum "candidate compression
feature" is characterized using its NATURAL bounding local minima in the
same gradient series (a data-determined extent, not a magnitude threshold).

A snapshot-relative display filter (peak magnitude > 0.5% of that snapshot's
own maximum gradient magnitude) is applied ONLY to suppress printing pure
round-off noise from the flat freestream region -- this is an analysis
display choice, explicitly NOT part of any final detection methodology.
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


def local_maxima(arr):
    """Indices i (1<=i<=len-2) where arr[i] >= both neighbors and > at
    least one neighbor (strict local maxima of the gradient-magnitude
    series). No magnitude threshold applied here."""
    idx = []
    for i in range(1, len(arr) - 1):
        if arr[i] >= arr[i - 1] and arr[i] >= arr[i + 1] and \
           (arr[i] > arr[i - 1] or arr[i] > arr[i + 1]):
            idx.append(i)
    return idx


def bounding_minima(arr, peak_i):
    """Walk outward from peak_i to the nearest local minima on each side
    (natural, data-determined extent of the compression feature)."""
    left = peak_i
    while left > 0 and arr[left - 1] <= arr[left]:
        left -= 1
    right = peak_i
    while right < len(arr) - 1 and arr[right + 1] <= arr[right]:
        right += 1
    return left, right


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
    ap.add_argument("--display-rel-floor", type=float, default=0.005,
                     help="Snapshot-relative display filter only "
                          "(fraction of that snapshot's max |grad|). "
                          "NOT a detection threshold.")
    args = ap.parse_args()

    case = Path(args.case_dir)
    poly = case / "constant" / "polyMesh"
    gamma = args.gamma
    R = args.R

    print("Loading mesh geometry...")
    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_labels(poly / "owner")
    neighbour = parse_labels(poly / "neighbour")
    centroids = cell_centroids_corrected(points, faces, owner, neighbour)
    print(f"  n_cells = {len(centroids)}")

    axis_mask = (np.hypot(centroids[:, 1], centroids[:, 2]) < args.axis_r_tol) & \
                (centroids[:, 0] > args.axis_x_min) & (centroids[:, 0] < args.axis_x_max)
    axis_idx_all = np.where(axis_mask)[0]
    axis_x_all = centroids[axis_idx_all, 0]
    order = np.argsort(axis_x_all)
    axis_idx_all = axis_idx_all[order]
    axis_x_all = axis_x_all[order]
    print(f"Axis cells: {len(axis_idx_all)}, x range "
          f"[{axis_x_all.min():.5f}, {axis_x_all.max():.5f}]\n")

    p2_p1_th, rho2_rho1_th, T2_T1_th = normal_shock_ratios(args.Minf, gamma)
    print(f"Theoretical normal-shock ratios (M_inf={args.Minf}, gamma={gamma}): "
          f"p2/p1={p2_p1_th:.3f}  rho2/rho1={rho2_rho1_th:.3f}  T2/T1={T2_T1_th:.3f}\n")

    tdirs = []
    for d in case.iterdir():
        if d.is_dir() and re.match(r'^[\d.eE+-]+$', d.name) and d.name != '0':
            try:
                t = float(d.name)
                tdirs.append((t, d))
            except ValueError:
                continue
    tdirs.sort(key=lambda x: x[0])

    if not tdirs:
        print("No snapshot directories found.")
        return

    rep_times = {tdirs[0][0]: "EARLY",
                 tdirs[len(tdirs) // 2][0]: "MIDDLE",
                 tdirs[-1][0]: "LATE"}

    all_snapshot_features = []  # (t, list of feature dicts) for compact summary

    for t, d in tdirs:
        p = parse_scalar_field(d / "p")
        rho = parse_scalar_field(d / "rho")
        T = parse_scalar_field(d / "T")
        U = parse_vector_field(d / "U")
        if p is None or rho is None or T is None or U is None:
            continue
        n = min(len(p), len(rho), len(T), len(U), len(centroids))
        valid = axis_idx_all < n
        aidx = axis_idx_all[valid]
        ax = axis_x_all[valid]
        if len(aidx) < 5:
            continue

        p_ax = p[aidx]
        rho_ax = rho[aidx]
        T_ax = T[aidx]
        Ux_ax = U[aidx, 0]
        a_ax = np.sqrt(gamma * R * np.maximum(T_ax, 1.0))
        M_ax = np.abs(Ux_ax) / a_ax

        dx = np.diff(ax)
        dpdx = np.abs(np.diff(p_ax) / dx)
        drhodx = np.abs(np.diff(rho_ax) / dx)
        # place gradient magnitude at midpoint index; for local-maxima
        # search treat as its own series
        xg = 0.5 * (ax[:-1] + ax[1:])

        peaks_p = local_maxima(dpdx)
        peaks_rho = local_maxima(drhodx)

        dpdx_max = dpdx.max() if len(dpdx) else 0.0
        drhodx_max = drhodx.max() if len(drhodx) else 0.0

        def characterize(peak_list, grad_arr, max_val, tag):
            feats = []
            for pi in peak_list:
                if max_val <= 0 or grad_arr[pi] < args.display_rel_floor * max_val:
                    continue
                li, ri = bounding_minima(grad_arr, pi)
                # map gradient-array indices back to axis-cell indices
                # grad_arr[i] corresponds to interval (ax[i], ax[i+1])
                left_cell = li
                right_cell = ri + 1
                x_peak = xg[pi]
                dist_from_nose = -x_peak
                dp_feat = p_ax[right_cell] - p_ax[left_cell]
                drho_feat = rho_ax[right_cell] - rho_ax[left_cell]
                dT_feat = T_ax[right_cell] - T_ax[left_cell]
                dM_feat = M_ax[right_cell] - M_ax[left_cell]
                n_cells = right_cell - left_cell
                feats.append(dict(tag=tag, x_peak=x_peak, dist=dist_from_nose,
                                   left_cell=left_cell, right_cell=right_cell,
                                   p_left=p_ax[left_cell], p_right=p_ax[right_cell],
                                   rho_left=rho_ax[left_cell], rho_right=rho_ax[right_cell],
                                   T_left=T_ax[left_cell], T_right=T_ax[right_cell],
                                   M_left=M_ax[left_cell], M_right=M_ax[right_cell],
                                   dp=dp_feat, drho=drho_feat, dT=dT_feat, dM=dM_feat,
                                   n_cells=n_cells, grad_val=grad_arr[pi]))
            return feats

        feats_p = characterize(peaks_p, dpdx, dpdx_max, "dp/dx")
        feats_rho = characterize(peaks_rho, drhodx, drhodx_max, "drho/dx")
        all_feats = sorted(feats_p + feats_rho, key=lambda f: f["x_peak"])
        all_snapshot_features.append((t, all_feats))

        if t in rep_times:
            label = rep_times[t]
            print("=" * 150)
            print(f"REPRESENTATIVE SNAPSHOT [{label}]  t = {t:.5e} s")
            print("=" * 150)
            print("Full axis profile (upstream -> body):")
            print(f"{'x':>10} {'p':>10} {'rho':>12} {'T':>9} {'M':>7}")
            for xi, pi, ri, Ti, Mi in zip(ax, p_ax, rho_ax, T_ax, M_ax):
                print(f"{xi:>10.5f} {pi:>10.3f} {ri:>12.5e} {Ti:>9.2f} {Mi:>7.3f}")
            print()
            print("Candidate compression features (sorted upstream -> body):")
            hdr = (f"{'src':>7} {'x_peak':>9} {'dist_nose':>10} {'cells':>6} "
                   f"{'p_L':>9} {'p_R':>10} {'rho_L':>11} {'rho_R':>11} "
                   f"{'T_L':>8} {'T_R':>8} {'M_L':>6} {'M_R':>6} {'gradval':>10}")
            print(hdr)
            for f in all_feats:
                print(f"{f['tag']:>7} {f['x_peak']:>9.5f} {f['dist']:>10.5f} "
                      f"{f['n_cells']:>6d} {f['p_left']:>9.2f} {f['p_right']:>10.2f} "
                      f"{f['rho_left']:>11.4e} {f['rho_right']:>11.4e} "
                      f"{f['T_left']:>8.2f} {f['T_right']:>8.2f} "
                      f"{f['M_left']:>6.2f} {f['M_right']:>6.2f} {f['grad_val']:>10.3e}")
            print()

    print("=" * 150)
    print("COMPACT CANDIDATE-FEATURE SUMMARY, ALL SNAPSHOTS "
          "(dp/dx-based candidates only, sorted upstream -> body per snapshot)")
    print("=" * 150)
    print(f"{'time':>12} {'src':>7} {'x_peak':>9} {'dist_nose':>10} {'cells':>6} "
          f"{'dp':>10} {'drho':>11} {'dT':>8} {'dM':>6}")
    for t, feats in all_snapshot_features:
        dp_feats = [f for f in feats if f["tag"] == "dp/dx"]
        for f in dp_feats:
            print(f"{t:>12.4e} {f['tag']:>7} {f['x_peak']:>9.5f} {f['dist']:>10.5f} "
                  f"{f['n_cells']:>6d} {f['dp']:>10.3f} {f['drho']:>11.4e} "
                  f"{f['dT']:>8.2f} {f['dM']:>6.3f}")


if __name__ == "__main__":
    main()
