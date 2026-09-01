"""
Project 01 - Hypersonic Blunt Body
Visual verification of Gmsh domain outline (Stage 2 mesh check).

Loads a .msh file, plots all mesh edges/triangles as a wireframe so the
domain shape, body profile, and farfield proportions can be visually
confirmed without relying on the Gmsh GUI (WSLg has shown some flakiness).
"""

import sys
import meshio
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def plot_mesh_outline(msh_path: Path, out_png: Path, title: str):
    mesh = meshio.read(str(msh_path))

    points = mesh.points[:, :2]  # x, y (r) only, drop z

    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot triangle mesh as a light wireframe (if present)
    for cell_block in mesh.cells:
        if cell_block.type == "triangle":
            tris = cell_block.data
            for tri in tris:
                pts = points[tri]
                pts_closed = np.vstack([pts, pts[0]])
                ax.plot(pts_closed[:, 0], pts_closed[:, 1], color="lightgray", lw=0.3, zorder=1)

    # Plot boundary lines distinctly (if present as line cells)
    for cell_block in mesh.cells:
        if cell_block.type == "line":
            lines = cell_block.data
            for ln in lines:
                pts = points[ln]
                ax.plot(pts[:, 0], pts[:, 1], color="tab:blue", lw=1.5, zorder=2)

    ax.set_xlabel("x [m] (axial)")
    ax.set_ylabel("r [m] (radial)")
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_png}")
    print(f"  Points: {len(points)}")
    for cb in mesh.cells:
        print(f"  {cb.type}: {len(cb.data)}")


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[3]
    proj = repo_root / "01-hypersonic-blunt-body"

    msh_path = proj / "mesh" / "baseline_domain_outline.msh"
    out_png = proj / "results" / "figures" / "domain_outline_check.png"

    plot_mesh_outline(msh_path, out_png, "Domain Outline Check (Stage 2, coarse, no BL)")
