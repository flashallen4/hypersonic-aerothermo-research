"""
Project 01 - Hypersonic Blunt Body
Zoomed visualization of near-wall mesh structure (boundary layer check).
"""

import meshio
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    repo_root = Path(__file__).resolve().parents[3]
    proj = repo_root / "01-hypersonic-blunt-body"
    msh_path = proj / "mesh" / "baseline_boundary_layer_check.msh"
    mesh = meshio.read(str(msh_path))
    points = mesh.points[:, :2]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    zoom_windows = [
        ("Stagnation region zoom (x: -0.002 to 0.008 m)", (-0.002, 0.008), (-0.002, 0.015)),
        ("Tangent point region zoom (x: 0.03 to 0.05 m)", (0.03, 0.05), (0.03, 0.06)),
    ]

    for ax, (title, xlim, ylim) in zip(axes, zoom_windows):
        for cell_block in mesh.cells:
            if cell_block.type in ("triangle", "quad"):
                for cell in cell_block.data:
                    pts = points[cell]
                    pts_closed = np.vstack([pts, pts[0]])
                    ax.plot(pts_closed[:, 0], pts_closed[:, 1], color="steelblue", lw=0.4)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("r [m]")

    out_png = proj / "results" / "figures" / "boundary_layer_zoom_check.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    print(f"Saved: {out_png}")

if __name__ == "__main__":
    main()
