import pyvista as pv
from dolfinx import plot

def create_live_plot(U0, c, dofsU):
    topology, cell_types, x = plot.vtk_mesh(U0)
    grid = pv.UnstructuredGrid(topology, cell_types, x)
    grid.point_data["u"] = c.x.array[dofsU].real
    grid.set_active_scalars("u")
    p = pv.Plotter(title="Attractiveness Field", auto_update=True)
    p.add_mesh(grid, opacity=0.8, cmap="viridis")
    p.view_xy(True)
    return p, grid

