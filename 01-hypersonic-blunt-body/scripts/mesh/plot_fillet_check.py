"""
Project 01 - Hypersonic Blunt Body
Visual verification of the cone/base fillet geometry: full-body view plus
a zoomed view of the fillet region, checking smooth tangent continuity.
"""

import meshio
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    repo_root = Path(__file__).resolve().parents[3]
    proj = repo_root / "01-hypersonic-blunt-body"
    msh_path = proj / "mesh" / "fillet_geometry_check.msh"
    mesh = meshio.read(str(msh_path))
    points = mesh.points[:, :2]

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    for cell_block in mesh.cells:
        if cell_block.type == "line":
            for ln in cell_block.data:
                pts = points[ln]
                axes[0].plot(pts[:, 0], pts[:, 1], color="tab:blue", lw=1.5)
                axes[1].plot(pts[:, 0], pts[:, 1], color="tab:blue", lw=1.5)

    axes[0].set_title("Full domain boundary")
    axes[0].set_xlabel("x [m]"); axes[0].set_ylabel("r [m]")
    axes[0].set_aspect("equal")
    axes[0].grid(alpha=0.3)

    L = 0.416622
    R_b = 0.15
    axes[1].set_xlim(L - 0.015, L + 0.003)
    axes[1].set_ylim(R_b - 0.015, R_b + 0.003)
    axes[1].set_title("Zoomed: cone/base fillet region")
    axes[1].set_xlabel("x [m]"); axes[1].set_ylabel("r [m]")
    axes[1].set_aspect("equal")
    axes[1].grid(alpha=0.3)

    out_png = proj / "results" / "figures" / "fillet_geometry_check.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    print(f"Saved: {out_png}")

if __name__ == "__main__":
    main()
