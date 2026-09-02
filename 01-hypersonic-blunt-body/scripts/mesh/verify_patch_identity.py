"""
Project 01 - Hypersonic Blunt Body
Directly verify OpenFOAM patch identity by parsing polyMesh points/faces
and computing face centroids for named patches, comparing against expected
geometric extents from CASE_SPECIFICATION.md.
"""

import re
from pathlib import Path
import numpy as np

def parse_foam_list(path, is_face_list=False):
    """Parse an OpenFOAM ASCII points or faces list file."""
    text = Path(path).read_text()
    # Strip comments
    text = re.sub(r'//.*', '', text)
    text = re.sub(r'/\*.*?\*/', '', text, flags=re.DOTALL)
    # Find the count and the parenthesized list
    lines = text.split('\n')
    # Find first standalone integer line (the count), then the following (...)
    joined = '\n'.join(lines)
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
    assert len(points) == count, f"{len(points)} != {count}"
    return points

def parse_faces(path):
    count, body = parse_foam_list(path)
    # Each face: N(p0 p1 p2 ...)
    faces = re.findall(r'\d+\(([^)]*)\)', body)
    face_list = [[int(x) for x in f.split()] for f in faces]
    assert len(face_list) == count, f"{len(face_list)} != {count}"
    return face_list

def main():
    repo_root = Path(__file__).resolve().parents[3]
    case = repo_root / "01-hypersonic-blunt-body" / "openfoam" / "wedge_test_case_v2"
    poly = case / "constant" / "polyMesh"

    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")

    boundary_text = (poly / "boundary").read_text()

    patches = {}
    for m in re.finditer(r'(\w+)\s*\{\s*type\s+\w+;\s*(?:physicalType\s+\w+;\s*)?nFaces\s+(\d+);\s*startFace\s+(\d+);', boundary_text):
        name, nfaces, startface = m.group(1), int(m.group(2)), int(m.group(3))
        patches[name] = (startface, nfaces)

    print(f"{'Patch':20s} {'nFaces':>7s} {'x_min':>10s} {'x_max':>10s} {'r_min':>10s} {'r_max':>10s}")
    for name, (start, n) in patches.items():
        centroids = []
        for fi in range(start, start + n):
            face_pts = points[faces[fi]]
            centroid = face_pts.mean(axis=0)
            centroids.append(centroid)
        centroids = np.array(centroids)
        x = centroids[:, 0]
        r = np.sqrt(centroids[:, 1]**2 + centroids[:, 2]**2)  # radial distance from axis
        print(f"{name:20s} {n:7d} {x.min():10.5f} {x.max():10.5f} {r.min():10.5f} {r.max():10.5f}")

    print("\nExpected geometric extents (baseline case, from CASE_SPECIFICATION.md):")
    print("  wall_nose:  x in [0, 0.03706],       r in [0, 0.04830]  (spherical cap)")
    print("  wall_cone:  x in [0.03706, 0.41662], r in [0.04830, 0.15]  (straight cone)")
    print("  wall_base:  x = 0.41662 (constant),  r in [0, 0.15]")

if __name__ == "__main__":
    main()
