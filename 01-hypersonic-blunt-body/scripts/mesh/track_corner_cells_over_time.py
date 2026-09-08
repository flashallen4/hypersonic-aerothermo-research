#!/usr/bin/env python3
"""
Track p/T/rho/U at specific cell IDs across all timestep directories in a
case, to assess whether the local flow state is converging toward
quasi-steady behavior or still actively evolving. Read-only.

Usage:
    python3 track_corner_cells_over_time.py <case_dir> --cell-ids 4632 4633 4634 4635 4636 4637
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--cell-ids", type=int, nargs="+", required=True)
    args = ap.parse_args()

    case = Path(args.case_dir)
    tdirs = []
    for d in case.iterdir():
        if d.is_dir() and re.match(r'^[\d.eE+-]+$', d.name):
            try:
                t = float(d.name)
                tdirs.append((t, d))
            except ValueError:
                continue
    tdirs.sort(key=lambda x: x[0])

    print(f"Found {len(tdirs)} timestep directories")
    print(f"Tracking cells: {args.cell_ids}\n")

    header = f"{'time':>14}"
    for cid in args.cell_ids:
        header += f" | {'p['+str(cid)+']':>12} {'T['+str(cid)+']':>10} {'rho['+str(cid)+']':>12} {'Ux['+str(cid)+']':>10}"
    print(header)
    print("-" * len(header))

    rows = []
    for t, d in tdirs:
        p, p_uni = parse_scalar_field(d / "p")
        T, T_uni = parse_scalar_field(d / "T")
        rho, rho_uni = parse_scalar_field(d / "rho")
        U = parse_vector_field(d / "U")

        line = f"{t:>14.6e}"
        row = {"time": t}
        for cid in args.cell_ids:
            pv = p[cid] if p is not None and cid < len(p) else (p_uni if p_uni is not None else float('nan'))
            Tv = T[cid] if T is not None and cid < len(T) else (T_uni if T_uni is not None else float('nan'))
            rv = rho[cid] if rho is not None and cid < len(rho) else (rho_uni if rho_uni is not None else float('nan'))
            uv = U[cid, 0] if U is not None and cid < len(U) else float('nan')
            line += f" | {pv:>12.4f} {Tv:>10.3f} {rv:>12.6e} {uv:>10.3f}"
            row[f"p_{cid}"] = pv
            row[f"T_{cid}"] = Tv
            row[f"rho_{cid}"] = rv
            row[f"Ux_{cid}"] = uv
        print(line)
        rows.append(row)

    print("\n--- Summary: min/max range across all snapshots per cell ---")
    for cid in args.cell_ids:
        ps = [r[f"p_{cid}"] for r in rows if not np.isnan(r[f"p_{cid}"])]
        Ts = [r[f"T_{cid}"] for r in rows if not np.isnan(r[f"T_{cid}"])]
        if ps:
            print(f"Cell {cid}: p range [{min(ps):.3f}, {max(ps):.3f}]  "
                  f"T range [{min(Ts):.3f}, {max(Ts):.3f}]")


if __name__ == "__main__":
    main()
