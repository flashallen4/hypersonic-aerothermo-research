#!/usr/bin/env python3
"""
Constructs 0/p, 0/T, 0/U for the isolated corner sub-case by mapping the
parent full-domain pre-crash snapshot onto:
  - the sub-mesh's internal field (nearest-neighbor, per sub-mesh cell)
  - cut_inlet_upstream boundary (nearest-neighbor, per face) -> fixedValue
  - cut_outer_radial boundary (nearest-neighbor, per face)   -> fixedValue
  - cut_outlet_downstream -> zeroGradient (no mapping needed)
  - wall_base / wall_cone -> same BCs as parent case (noSlip, T=300, p zeroGradient)
  - frontWedge / backWedge -> wedge

Also reports mapping-quality diagnostics (distances, duplicates, field
min/max) before any file is written, so the mapping can be inspected.

Does NOT modify the parent case. Does NOT run the solver. Does NOT change
any numerical/solver settings.
"""

import re
import argparse
from pathlib import Path
from collections import Counter
import numpy as np
from scipy.spatial import cKDTree


# ---------- shared parsers (validated earlier in session, reused) ----------

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


def cell_centroids_corrected(points, faces, owner, neighbour):
    """Corrected: accumulates each internal face's centre into BOTH the
    owner and neighbour cell (previous scripts only used owner)."""
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
            n_ = neighbour[fid]
            csum[n_] += fc
            ccount[n_] += 1
    return csum / np.maximum(ccount, 1)[:, None]


# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent-case", required=True)
    ap.add_argument("--parent-timestep", required=True)
    ap.add_argument("--sub-case", required=True)
    args = ap.parse_args()

    parent = Path(args.parent_case)
    sub = Path(args.sub_case)
    ptdir = parent / args.parent_timestep

    print(f"Parent case: {parent}, timestep: {args.parent_timestep}")
    print(f"Sub-case: {sub}")

    # --- Parent mesh + fields ---
    p_points = parse_points(parent / "constant/polyMesh/points")
    p_faces = parse_faces(parent / "constant/polyMesh/faces")
    p_owner = parse_labels(parent / "constant/polyMesh/owner")
    p_neighbour = parse_labels(parent / "constant/polyMesh/neighbour")

    print("Computing parent cell centroids (corrected owner+neighbour method)...")
    parent_centroids = cell_centroids_corrected(p_points, p_faces, p_owner, p_neighbour)
    print(f"  parent n_cells = {len(parent_centroids)}")

    T, T_uni = parse_scalar_field(ptdir / "T")
    p, p_uni = parse_scalar_field(ptdir / "p")
    rho, rho_uni = parse_scalar_field(ptdir / "rho")
    U = parse_vector_field(ptdir / "U")

    tree = cKDTree(parent_centroids)

    # --- Sub-mesh ---
    s_points = parse_points(sub / "constant/polyMesh/points")
    s_faces = parse_faces(sub / "constant/polyMesh/faces")
    s_owner = parse_labels(sub / "constant/polyMesh/owner")
    s_neighbour = parse_labels(sub / "constant/polyMesh/neighbour")
    s_boundary = parse_boundary(sub / "constant/polyMesh/boundary")

    print("Computing sub-mesh cell centroids...")
    sub_centroids = cell_centroids_corrected(s_points, s_faces, s_owner, s_neighbour)
    n_sub_cells = len(sub_centroids)
    print(f"  sub-mesh n_cells = {n_sub_cells}")

    def patch_info(name):
        for n, t, nf, sf in s_boundary:
            if n == name:
                return nf, sf
        raise KeyError(name)

    def face_centroids_for_patch(name):
        nf, sf = patch_info(name)
        cens = []
        for fid in range(sf, sf + nf):
            cens.append(s_points[s_faces[fid]].mean(axis=0))
        return np.array(cens)

    # --- Map internal field (every sub-mesh cell) ---
    print("\nMapping internal field (sub-mesh cells -> nearest parent cell)...")
    int_dist, int_idx = tree.query(sub_centroids)
    print(f"  internal mapping distance: min={int_dist.min():.3e} "
          f"mean={int_dist.mean():.3e} max={int_dist.max():.3e}")
    dup_counts = Counter(int_idx.tolist())
    n_dup = sum(1 for v in dup_counts.values() if v > 1)
    print(f"  unique parent cells used: {len(dup_counts)} / {n_sub_cells} sub-cells "
          f"({n_dup} parent cells reused by >1 sub-cell)")

    T_int = T[int_idx]
    p_int = p[int_idx]
    rho_int = rho[int_idx]
    U_int = U[int_idx]

    Cv = 718.0
    Tstd = 298.15
    e_int = Cv * (T_int - Tstd)

    # --- Map cut_inlet_upstream ---
    print("\nMapping cut_inlet_upstream faces...")
    inlet_cens = face_centroids_for_patch("cut_inlet_upstream")
    inlet_dist, inlet_idx = tree.query(inlet_cens)
    print(f"  n faces = {len(inlet_cens)}")
    print(f"  mapping distance: min={inlet_dist.min():.3e} mean={inlet_dist.mean():.3e} "
          f"max={inlet_dist.max():.3e}")
    dup = Counter(inlet_idx.tolist())
    print(f"  unique parent cells used: {len(dup)} / {len(inlet_cens)} "
          f"({sum(1 for v in dup.values() if v > 1)} reused)")
    T_inlet = T[inlet_idx]; p_inlet = p[inlet_idx]
    rho_inlet = rho[inlet_idx]; U_inlet = U[inlet_idx]
    e_inlet = Cv * (T_inlet - Tstd)
    print(f"  T: min={T_inlet.min():.3f} max={T_inlet.max():.3f}")
    print(f"  p: min={p_inlet.min():.3f} max={p_inlet.max():.3f}")
    print(f"  rho: min={rho_inlet.min():.6e} max={rho_inlet.max():.6e}")
    print(f"  Ux: min={U_inlet[:,0].min():.3f} max={U_inlet[:,0].max():.3f}")
    print(f"  e(diagnostic, Cv*(T-Tstd)): min={e_inlet.min():.3f} max={e_inlet.max():.3f}")

    # --- Map cut_outer_radial ---
    print("\nMapping cut_outer_radial faces...")
    outer_cens = face_centroids_for_patch("cut_outer_radial")
    outer_dist, outer_idx = tree.query(outer_cens)
    print(f"  n faces = {len(outer_cens)}")
    print(f"  mapping distance: min={outer_dist.min():.3e} mean={outer_dist.mean():.3e} "
          f"max={outer_dist.max():.3e}")
    dup = Counter(outer_idx.tolist())
    print(f"  unique parent cells used: {len(dup)} / {len(outer_cens)} "
          f"({sum(1 for v in dup.values() if v > 1)} reused)")
    T_outer = T[outer_idx]; p_outer = p[outer_idx]
    rho_outer = rho[outer_idx]; U_outer = U[outer_idx]
    e_outer = Cv * (T_outer - Tstd)
    print(f"  T: min={T_outer.min():.3f} max={T_outer.max():.3f}")
    print(f"  p: min={p_outer.min():.3f} max={p_outer.max():.3f}")
    print(f"  rho: min={rho_outer.min():.6e} max={rho_outer.max():.6e}")
    print(f"  Ux: min={U_outer[:,0].min():.3f} max={U_outer[:,0].max():.3f}")
    print(f"  e(diagnostic, Cv*(T-Tstd)): min={e_outer.min():.3f} max={e_outer.max():.3f}")

    # --- Write field files ---
    def scalar_field_text(name, dims, internal_vals, patch_blocks):
        lines = []
        lines.append('/*--------------------------------*- C++ -*----------------------------------*\\')
        lines.append('\\*---------------------------------------------------------------------------*/')
        lines.append('FoamFile\n{\n    version 2.0;\n    format ascii;\n    class volScalarField;\n'
                      f'    object {name};\n}}')
        lines.append(f'dimensions {dims};')
        lines.append(f'internalField nonuniform List<scalar>')
        lines.append(str(len(internal_vals)))
        lines.append('(')
        lines.extend(f'{v:.6g}' for v in internal_vals)
        lines.append(');')
        lines.append('boundaryField')
        lines.append('{')
        lines.extend(patch_blocks)
        lines.append('}')
        return "\n".join(lines) + "\n"

    def vector_field_text(name, dims, internal_vals, patch_blocks):
        lines = []
        lines.append('/*--------------------------------*- C++ -*----------------------------------*\\')
        lines.append('\\*---------------------------------------------------------------------------*/')
        lines.append('FoamFile\n{\n    version 2.0;\n    format ascii;\n    class volVectorField;\n'
                      f'    object {name};\n}}')
        lines.append(f'dimensions {dims};')
        lines.append(f'internalField nonuniform List<vector>')
        lines.append(str(len(internal_vals)))
        lines.append('(')
        lines.extend(f'({v[0]:.6g} {v[1]:.6g} {v[2]:.6g})' for v in internal_vals)
        lines.append(');')
        lines.append('boundaryField')
        lines.append('{')
        lines.extend(patch_blocks)
        lines.append('}')
        return "\n".join(lines) + "\n"

    def nonuniform_scalar_patch(name, ptype, vals):
        s = f"    {name}\n    {{\n        type {ptype};\n"
        if vals is not None:
            s += f"        value nonuniform List<scalar>\n        {len(vals)}\n        (\n"
            s += "\n".join(f"        {v:.6g}" for v in vals)
            s += "\n        );\n"
        s += "    }\n"
        return s

    def nonuniform_vector_patch(name, ptype, vals):
        s = f"    {name}\n    {{\n        type {ptype};\n"
        if vals is not None:
            s += f"        value nonuniform List<vector>\n        {len(vals)}\n        (\n"
            s += "\n".join(f"        ({v[0]:.6g} {v[1]:.6g} {v[2]:.6g})" for v in vals)
            s += "\n        );\n"
        s += "    }\n"
        return s

    def uniform_patch(name, ptype, value=None):
        if value is None:
            return f"    {name}\n    {{\n        type {ptype};\n    }}\n"
        return f"    {name}\n    {{\n        type {ptype};\n        value uniform {value};\n    }}\n"

    zero_dir = sub / "0"
    zero_dir.mkdir(exist_ok=True)

    # p
    p_patches = [
        nonuniform_scalar_patch("cut_inlet_upstream", "fixedValue", p_inlet),
        uniform_patch("cut_outlet_downstream", "zeroGradient"),
        nonuniform_scalar_patch("cut_outer_radial", "fixedValue", p_outer),
        uniform_patch("wall_base", "zeroGradient"),
        uniform_patch("wall_cone", "zeroGradient"),
        uniform_patch("frontWedge", "wedge"),
        uniform_patch("backWedge", "wedge"),
        uniform_patch("oldInternalFaces", "internal"),
    ]
    (zero_dir / "p").write_text(scalar_field_text("p", "[1 -1 -2 0 0 0 0]", p_int, p_patches))

    # T
    T_patches = [
        nonuniform_scalar_patch("cut_inlet_upstream", "fixedValue", T_inlet),
        uniform_patch("cut_outlet_downstream", "zeroGradient"),
        nonuniform_scalar_patch("cut_outer_radial", "fixedValue", T_outer),
        uniform_patch("wall_base", "fixedValue", 300),
        uniform_patch("wall_cone", "fixedValue", 300),
        uniform_patch("frontWedge", "wedge"),
        uniform_patch("backWedge", "wedge"),
        uniform_patch("oldInternalFaces", "internal"),
    ]
    (zero_dir / "T").write_text(scalar_field_text("T", "[0 0 0 1 0 0 0]", T_int, T_patches))

    # U
    U_patches = [
        nonuniform_vector_patch("cut_inlet_upstream", "fixedValue", U_inlet),
        uniform_patch("cut_outlet_downstream", "zeroGradient"),
        nonuniform_vector_patch("cut_outer_radial", "fixedValue", U_outer),
        uniform_patch("wall_base", "noSlip"),
        uniform_patch("wall_cone", "noSlip"),
        uniform_patch("frontWedge", "wedge"),
        uniform_patch("backWedge", "wedge"),
        uniform_patch("oldInternalFaces", "internal"),
    ]
    (zero_dir / "U").write_text(vector_field_text("U", "[0 1 -1 0 0 0 0]", U_int, U_patches))

    print(f"\nWrote {zero_dir}/p, {zero_dir}/T, {zero_dir}/U")
    print(f"internalField: {n_sub_cells} cells mapped")
    print(f"cut_inlet_upstream: {len(inlet_cens)} faces mapped")
    print(f"cut_outer_radial: {len(outer_cens)} faces mapped")
    print("cut_outlet_downstream: zeroGradient (no mapping)")


if __name__ == "__main__":
    main()
