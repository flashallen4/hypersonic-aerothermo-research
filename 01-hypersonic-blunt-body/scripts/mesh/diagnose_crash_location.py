"""
Project 01 - Hypersonic Blunt Body
Diagnose solver crash location by inspecting the last-written field state
for extreme/suspicious values (T, p, U magnitude), and reporting the
corresponding cell-centre coordinates.
"""

import re
from pathlib import Path
import numpy as np

def parse_scalar_field(path):
    text = Path(path).read_text()
    text_nc = re.sub(r'//.*', '', text)
    m = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\n?\s*(\d+)\s*\n\s*\(([^)]*)\)', text_nc, re.DOTALL)
    if m is None:
        # uniform field
        m2 = re.search(r'internalField\s+uniform\s+([\-\d.eE]+)', text_nc)
        return None
    n = int(m.group(1))
    vals = np.array([float(x) for x in m.group(2).split()])
    assert len(vals) == n
    return vals

def parse_vector_field(path):
    text = Path(path).read_text()
    text_nc = re.sub(r'//.*', '', text)
    m = re.search(r'internalField\s+nonuniform\s+List<vector>\s*\n?\s*(\d+)\s*\n\s*\(([^)]*)\)', text_nc, re.DOTALL)
    if m is None:
        return None
    n = int(m.group(1))
    tuples = re.findall(r'\(([^)]+)\)', m.group(2))
    vecs = np.array([[float(x) for x in t.split()] for t in tuples])
    assert len(vecs) == n
    return vecs

def parse_cell_centres(case_dir):
    """Approximate cell centres via mean of face centres from points+faces+owner+neighbour."""
    poly = case_dir / "constant" / "polyMesh"
    # Reuse the parser pattern from verify_patch_identity.py
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
        points = np.array([[float(x) for x in p.split()] for p in pts])
        return points

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

    cell_centres = cell_face_sum / np.maximum(cell_face_count, 1)[:, None]
    return cell_centres

def main():
    repo_root = Path(__file__).resolve().parents[3]
    case = repo_root / "01-hypersonic-blunt-body" / "openfoam" / "baseline_case"
    last_time = "1.41839e-08"
    tdir = case / last_time

    print(f"Loading fields from {tdir}")
    T = parse_scalar_field(tdir / "T")
    p = parse_scalar_field(tdir / "p")
    U = parse_vector_field(tdir / "U")

    print(f"T: n={None if T is None else len(T)}, min={T.min() if T is not None else 'NA'}, max={T.max() if T is not None else 'NA'}")
    print(f"p: n={None if p is None else len(p)}, min={p.min() if p is not None else 'NA'}, max={p.max() if p is not None else 'NA'}")
    if U is not None:
        Umag = np.linalg.norm(U, axis=1)
        print(f"|U|: min={Umag.min()}, max={Umag.max()}")

    print("\nComputing cell centres (this may take a moment)...")
    centres = parse_cell_centres(case)
    print(f"n_cells (from mesh) = {len(centres)}")

    if T is not None:
        idx_sorted = np.argsort(T)[::-1]  # highest T first
        print("\n--- Top 10 highest-T cells ---")
        for i in idx_sorted[:10]:
            x, y, z = centres[i]
            r = np.sqrt(y**2 + z**2)
            print(f"  cell {i}: T={T[i]:.2f} K, x={x:.5f}, r={r:.5f}")

        idx_sorted_low = np.argsort(T)
        print("\n--- Top 10 LOWEST-T cells (check for negative/nonphysical) ---")
        for i in idx_sorted_low[:10]:
            x, y, z = centres[i]
            r = np.sqrt(y**2 + z**2)
            print(f"  cell {i}: T={T[i]:.2f} K, x={x:.5f}, r={r:.5f}")

if __name__ == "__main__":
    main()
