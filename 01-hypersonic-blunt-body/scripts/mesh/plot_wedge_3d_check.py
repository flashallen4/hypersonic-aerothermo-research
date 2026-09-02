"""
Project 01 - Hypersonic Blunt Body
3D visualization of wedge extrusion test mesh (Stage 4a).

Loads via meshio directly (not pv.read/from_meshio) to avoid a known
meshio/pyvista cell_sets_to_data() indexing bug triggered by Gmsh's
default entity tagging on extruded meshes without explicit Physical Groups.
"""

import meshio
import numpy as np
import pyvista as pv
from pathlib import Path


def meshio_to_pyvista_unstructured(mesh: meshio.Mesh) -> pv.UnstructuredGrid:
    """Manually build a PyVista UnstructuredGrid from meshio cell data,
    ignoring cell_sets/physical-group metadata entirely."""
    points = mesh.points

    # VTK cell type codes for the element types we expect from a wedge extrusion
    vtk_type_map = {
        "tetra": 10,
        "hexahedron": 12,
        "wedge": 13,   # prism
        "pyramid": 14,
        "triangle": 5,
        "quad": 9,
        "line": 3,
        "vertex": 1,
    }
    n_pts_map = {
        "tetra": 4, "hexahedron": 8, "wedge": 6, "pyramid": 5,
        "triangle": 3, "quad": 4, "line": 2, "vertex": 1,
    }

    # Keep only 3D volumetric cells for the solid-volume view
    volumetric_types = {"tetra", "hexahedron", "wedge", "pyramid"}

    cells = []
    cell_types = []
    for cell_block in mesh.cells:
        if cell_block.type not in volumetric_types:
            continue
        n_per = n_pts_map[cell_block.type]
        vtk_type = vtk_type_map[cell_block.type]
        for conn in cell_block.data:
            cells.append(n_per)
            cells.extend(conn.tolist())
            cell_types.append(vtk_type)

    cells = np.array(cells, dtype=np.int64)
    cell_types = np.array(cell_types, dtype=np.uint8)

    grid = pv.UnstructuredGrid(cells, cell_types, points)
    return grid


def main():
    repo_root = Path(__file__).resolve().parents[3]
    proj = repo_root / "01-hypersonic-blunt-body"
    msh_path = proj / "mesh" / "wedge_extrude_test.msh"

    mesh = meshio.read(str(msh_path))
    print("Cell blocks in file:")
    for cb in mesh.cells:
        print(f"  {cb.type}: {len(cb.data)}")

    grid = meshio_to_pyvista_unstructured(mesh)
    print(f"\nPyVista grid: {grid.n_cells} volumetric cells, {grid.n_points} points")
    print(f"Bounds: {grid.bounds}")

    plotter = pv.Plotter(off_screen=True, window_size=[1400, 900])
    plotter.add_mesh(grid, show_edges=True, color="lightblue", opacity=1.0)
    plotter.camera_position = [
        (0.8, 0.5, 0.3),
        (0.1, 0.05, 0.0),
        (0, 0, 1),
    ]
    plotter.add_axes()

    out_png = proj / "results" / "figures" / "wedge_extrude_3d_check.png"
    plotter.screenshot(str(out_png))
    print(f"\nSaved: {out_png}")


if __name__ == "__main__":
    main()
