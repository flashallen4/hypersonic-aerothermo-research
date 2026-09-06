#!/usr/bin/env python3
"""
Generalized diagnostic tool (v2): fixes the vector-field parser bug from
diagnose_crash_location_generalized.py (naive [^)]* regex truncated at the
first vector's closing paren). Uses balanced-paren extraction instead.
Read-only: parses field files and mesh geometry, does not modify anything.
Accepts case directory and timestep as arguments.
"""

import re
import argparse
from pathlib import Path
import numpy as np


def find_balanced_block(text, start_idx):
    """Given text and the index of an opening '(', return the substring
    between it and its matching closing ')', plus the index just after
    the closing paren."""
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
    raise ValueError("Unbalanced parentheses")


def parse_scalar_field(path):
    text = Path(path).read_text()
    text_nc = re.sub(r'//.*', '', text)
    m = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\n?\s*(\d+)\s*\n', text_nc)
    if m is None:
        m2 = re.search(r'internalField\s+uniform\s+([\-\d.eE]+)', text_nc)
        if m2:
            return None, float(m2.group(1))
        return None, None
    n = int(m.group(1))
    # find the opening paren right after the count
    open_idx = text_nc.index('(', m.end() - 1)
    body, _ = find_balanced_block(text_nc, open_idx)
    vals = np.array([float(x) for x in body.split()])
    assert len(vals) == n, f"expected {n}, got {len(vals)}"
    return vals, None


def parse_vector_field(path):
    """CORRECTED: uses balanced-paren extraction for the outer list, then
    finds each individual (x y z) tuple within it."""
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
    assert len(vecs) == n, f"expected {n}, got {len(vecs)}"
    return vecs


def parse_cell_centres(case_dir):
    poly = case_dir / "constant" / "polyMesh"

    def parse_foam_list(path):
        text = Path(path).read_text()
        text = re.sub(r'//.*', '', text)
        text = re.sub(r'/\*.*?\*/', '', text, flags=re.DOTALL)
        joined = text
        m = re.search(r'\n\s*(\d+)\s*\n\s*\(', joined)
        count = int(m.group(1))
        start = m.end()
        depth = 1
        i = start
        while depth > 0:
            if joined[i] == '(':
                depth += 1
            elif joined[i] == ')':
                depth -= 1
            i += 1
        body = joined[start:i-1]
        return count, body

    def parse_points(path):
        count, body = parse_foam_list(path)
        pts = re.findall(r'\(([^)]+)\)', body)
        return np.array([[float(x) for x in p.split()] for p in pts])

    def parse_faces(path):
        count, body = parse_foam_list(path)
        faces = re.findall(r'\d+\(([^)]*)\)', body)
        return [[int(x) for x in f.split()] for f in faces]

    def parse_owner(path):
        count, body = parse_foam_list(path)
        return np.array([int(x) for x in body.split()])

    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_owner(poly / "owner")

    n_cells = owner.max() + 1
    cell_face_sum = np.zeros((n_cells, 3))
    cell_face_count = np.zeros(n_cells)
    for fi, face in enumerate(faces):
        if fi >= len(owner):
            break
        face_centre = points[face].mean(axis=0)
        c = owner[fi]
        cell_face_sum[c] += face_centre
        cell_face_count[c] += 1

    return cell_face_sum / np.maximum(cell_face_count, 1)[:, None]


def report_scalar(name, arr, uniform_val, centres, k=5):
    print(f"\n--- {name} ---")
    if arr is None and uniform_val is not None:
        print(f"  UNIFORM field, value={uniform_val} (no spatial variation)")
        return
    if arr is None:
        print(f"  Could not parse field.")
        return
    print(f"  n={len(arr)}, min={arr.min():.6g}, max={arr.max():.6g}, mean={arr.mean():.6g}")
    idx_hi = np.argsort(arr)[::-1][:k]
    idx_lo = np.argsort(arr)[:k]
    print(f"  Top {k} HIGHEST-value cells:")
    for i in idx_hi:
        x, y, z = centres[i]
        r = np.hypot(y, z)
        print(f"    cell {i}: value={arr[i]:.6g}  x={x:.6f} r={r:.6f}")
    print(f"  Top {k} LOWEST-value cells:")
    for i in idx_lo:
        x, y, z = centres[i]
        r = np.hypot(y, z)
        print(f"    cell {i}: value={arr[i]:.6g}  x={x:.6f} r={r:.6f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("timestep")
    args = ap.parse_args()

    case = Path(args.case_dir)
    tdir = case / args.timestep
    print(f"Loading fields from {tdir}")

    T, T_uni = parse_scalar_field(tdir / "T")
    p, p_uni = parse_scalar_field(tdir / "p")
    rho, rho_uni = (parse_scalar_field(tdir / "rho")
                    if (tdir / "rho").exists() else (None, None))
    U = parse_vector_field(tdir / "U") if (tdir / "U").exists() else None

    print("\nComputing cell centres...")
    centres = parse_cell_centres(case)
    print(f"n_cells (from mesh) = {len(centres)}")

    report_scalar("T [K]", T, T_uni, centres)
    report_scalar("p [Pa]", p, p_uni, centres)
    report_scalar("rho [kg/m^3]", rho, rho_uni, centres)

    if U is not None:
        Umag = np.linalg.norm(U, axis=1)
        Ux = U[:, 0]
        report_scalar("|U| [m/s]", Umag, None, centres)
        report_scalar("Ux [m/s]", Ux, None, centres)
    else:
        print("\n--- U ---\n  Not found or could not parse.")

    print("\n--- e (sensibleInternalEnergy) ---")
    if (tdir / "e").exists():
        e, e_uni = parse_scalar_field(tdir / "e")
        report_scalar("e [J/kg]", e, e_uni, centres)
    else:
        print("  NOT WRITTEN at this timestep - cannot inspect directly.")

    # Targeted region check: x~0.416, r~0.149
    print(f"\n{'='*70}")
    print(f"TARGETED REGION CHECK: x in [0.415,0.418], r in [0.147,0.151]")
    print(f"{'='*70}")
    mask = ((centres[:, 0] > 0.415) & (centres[:, 0] < 0.418) &
            (np.hypot(centres[:, 1], centres[:, 2]) > 0.147) &
            (np.hypot(centres[:, 1], centres[:, 2]) < 0.151))
    idxs = np.where(mask)[0]
    print(f"Cells in this region: {len(idxs)}")
    if len(idxs) > 0 and T is not None:
        print(f"  T in region: min={T[idxs].min():.4g} max={T[idxs].max():.4g}")
    if len(idxs) > 0 and p is not None:
        print(f"  p in region: min={p[idxs].min():.4g} max={p[idxs].max():.4g}")
    if len(idxs) > 0 and rho is not None:
        print(f"  rho in region: min={rho[idxs].min():.4g} max={rho[idxs].max():.4g}")
    if len(idxs) > 0 and U is not None:
        Ux_r = U[idxs, 0]
        print(f"  Ux in region: min={Ux_r.min():.4g} max={Ux_r.max():.4g}")


if __name__ == "__main__":
    main()
