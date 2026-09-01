"""
Project 01 - Hypersonic Blunt Body
Full baseline mesh verification (Stage 3b): domain-wide view + shock-zone
zoom + boundary-layer zoom, to confirm both refinement mechanisms coexist
correctly.
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
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    ax.set_aspect("equal")


def main():
    repo_root = Path(__file__).resolve().parents[3]
    proj = repo_root / "01-hypersonic-blunt-body"
    msh_path = proj / "mesh" / "baseline_mesh_full.msh"
    mesh = meshio.read(str(msh_path))
    points = mesh.points[:, :2]

    R_n = 0.05
    theta_c = np.radians(15.0)
    R_b = 0.15
    x_t = R_n * (1 - np.sin(theta_c))
    r_t = R_n * np.cos(theta_c)
    L = x_t + (R_b - r_t) / np.tan(theta_c)

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # Panel 1: full domain
    plot_edges(axes[0], mesh, points, lw=0.15)
    axes[0].set_title("Full domain")
    axes[0].set_xlabel("x [m]"); axes[0].set_ylabel("r [m]")

    # Panel 2: shock-zone scale view (body + surrounding refined region)
    plot_edges(axes[1], mesh, points, xlim=(-0.2, 0.6), ylim=(-0.02, 0.4), lw=0.2)
    axes[1].set_title("Shock-layer refinement zone view")
    axes[1].set_xlabel("x [m]"); axes[1].set_ylabel("r [m]")

    # Panel 3: near-wall zoom (confirm BL field still intact)
    plot_edges(axes[2], mesh, points, xlim=(-0.002, 0.01), ylim=(-0.002, 0.02), lw=0.3)
    axes[2].set_title("Stagnation-region zoom (BL check)")
    axes[2].set_xlabel("x [m]"); axes[2].set_ylabel("r [m]")

    out_png = proj / "results" / "figures" / "full_mesh_check.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    print(f"Saved: {out_png}")
    print(f"Total nodes: {len(points)}")
    for cb in mesh.cells:
        print(f"  {cb.type}: {len(cb.data)}")


if __name__ == "__main__":
    main()
