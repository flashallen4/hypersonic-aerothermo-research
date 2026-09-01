"""
Project 01 - Hypersonic Blunt Body
Rigorous check of near-wall mesh spacing against boundary-layer sizing targets.

Method: build a dense analytic reference curve (matching blunted_cone_profile.py
exactly), then for every mesh node compute the perpendicular distance to the
nearest reference curve point. Bucket by arc-length position along the wall,
and report the smallest nonzero distance in each bucket as the achieved
first-layer wall-normal height at that station. This avoids the flaw of
"nearest node to a hand-picked approximate point," which conflates tangential
wall-vertex spacing with normal boundary-layer spacing.
"""

import meshio
import numpy as np
from pathlib import Path

TARGET_FIRST_CELL = 5.607e-6  # m

R_n = 0.05
theta_c = np.radians(15.0)
R_b = 0.15
x_t = R_n * (1 - np.sin(theta_c))
r_t = R_n * np.cos(theta_c)
L = x_t + (R_b - r_t) / np.tan(theta_c)


def analytic_wall_curve(n_nose=4000, n_cone=4000):
    x_nose = np.linspace(0.0, x_t, n_nose)
    r_nose = np.sqrt(np.clip(R_n**2 - (x_nose - R_n) ** 2, 0.0, None))
    x_cone = np.linspace(x_t, L, n_cone)[1:]
    r_cone = r_t + (x_cone - x_t) * np.tan(theta_c)
    x = np.concatenate([x_nose, x_cone])
    r = np.concatenate([r_nose, r_cone])
    return x, r


def main():
    repo_root = Path(__file__).resolve().parents[3]
    proj = repo_root / "01-hypersonic-blunt-body"
    msh_path = proj / "mesh" / "baseline_boundary_layer_check.msh"

    mesh = meshio.read(str(msh_path))
    points = mesh.points[:, :2]

    x_ref, r_ref = analytic_wall_curve()
    ref_pts = np.column_stack([x_ref, r_ref])

    # Only consider mesh nodes reasonably close to the wall (within 0.01 m)
    # to keep the nearest-reference-point search cheap and relevant.
    near_wall_mask = points[:, 1] < 0.20  # generous radial cutoff near the body
    candidate_idx = np.where(near_wall_mask)[0]
    candidates = points[candidate_idx]

    # For each candidate mesh node, find nearest reference curve point
    # (brute-force in chunks; candidate count should be modest near the wall)
    from scipy.spatial import cKDTree
    tree = cKDTree(ref_pts)
    dist_to_wall, nearest_ref_idx = tree.query(candidates)

    # Stations to inspect, by reference-curve index (approx arc-length position)
    stations = {
        "near stagnation (x~0.002)": np.argmin(np.abs(x_ref - 0.002)),
        "mid-nose (x~0.025)": np.argmin(np.abs(x_ref - 0.025)),
        "tangent region (x~x_t)": np.argmin(np.abs(x_ref - x_t)),
        "mid-cone (x~0.2)": np.argmin(np.abs(x_ref - 0.2)),
        "near base (x~0.40)": np.argmin(np.abs(x_ref - 0.40)),
    }

    print(f"Target first-cell height: {TARGET_FIRST_CELL:.4e} m ({TARGET_FIRST_CELL*1e6:.3f} um)\n")

    window = 60  # reference-point index window (+/-) to bucket nearby mesh nodes
    for label, ref_idx in stations.items():
        in_bucket = np.abs(nearest_ref_idx - ref_idx) <= window
        bucket_dists = dist_to_wall[in_bucket]
        bucket_dists = np.sort(bucket_dists)
        # smallest nonzero distance = achieved first-layer height at this station
        nonzero = bucket_dists[bucket_dists > 1e-9]
        if len(nonzero) == 0:
            print(f"{label:28s}: no nonzero-distance nodes found in bucket (n={len(bucket_dists)})")
            continue
        achieved = nonzero[0]
        ratio = achieved / TARGET_FIRST_CELL
        flag = "  <-- WAY OFF" if (ratio > 5 or ratio < 0.2) else ""
        print(f"{label:28s}: achieved = {achieved:.4e} m ({achieved*1e6:8.3f} um)  ratio={ratio:6.2f}{flag}  [n={len(bucket_dists)}]")


if __name__ == "__main__":
    main()
