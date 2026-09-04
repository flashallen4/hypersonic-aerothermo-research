"""
Project 01 - Hypersonic Blunt Body
Visual check of the adopted Case 4 baseline mesh (sharp corner, extended
EdgesList): full domain, cone/base corner zoom (confirm defect-free), and
the newly-located skew face region (x~0.55, r~0.005, wake/axis area).
"""

import meshio
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def plot_edges(ax, mesh, points, xlim=None, ylim=None, lw=0.3, color="steelblue"):
    for cell_block in mesh.cells:
        if cell_block.type in ("triangle", "quad"):
            for cell in cell_block.data:
                pts = points[cell]
                pts_closed = np.vstack([pts, pts[0]])
                ax.plot(pts_closed[:, 0], pts_closed[:, 1], color=color, lw=lw)
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    ax.set_aspect("equal")

def main():
    repo_root = Path(__file__).resolve().parents[3]
    proj = repo_root / "01-hypersonic-blunt-body"
    msh_path = proj / "mesh" / "baseline_mesh_wedge.msh"
    mesh = meshio.read(str(msh_path))
    points = mesh.points[:, :2]

    fig, axes = plt.subplots(1, 3, figsize=(21, 6))

    plot_edges(axes[0], mesh, points, lw=0.15)
    axes[0].set_title("Full domain (Case 4: sharp corner, extended EdgesList)")
    axes[0].set_xlabel("x [m]"); axes[0].set_ylabel("r [m]")

    plot_edges(axes[1], mesh, points, xlim=(0.40, 0.43), ylim=(0.13, 0.16), lw=0.3)
    axes[1].set_title("Cone/base corner zoom (should be defect-free)")
    axes[1].set_xlabel("x [m]"); axes[1].set_ylabel("r [m]")

    plot_edges(axes[2], mesh, points, xlim=(0.50, 0.60), ylim=(-0.01, 0.02), lw=0.3)
    axes[2].set_title("Skew face region zoom (x~0.55, r~0.005, wake/axis)")
    axes[2].set_xlabel("x [m]"); axes[2].set_ylabel("r [m]")

    out_png = proj / "results" / "figures" / "baseline_case4_check.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    print(f"Saved: {out_png}")

if __name__ == "__main__":
    main()
