import gmsh
import numpy as np
from mpi4py import MPI
from dolfinx import mesh
from dolfinx.io import gmshio
from dolfinx.fem import functionspace
from basix.ufl import element, mixed_element




def create_domain(params):
    """Create rectangular mesh and mixed function spaces."""
    Lx,Ly = params["Lx"], params["Ly"]
    nx, ny = params["Nx"], params["Ny"]
    corner_1 = np.array([0, 0])
    corner_2 = np.array([Lx, Ly])
    msh = mesh.create_rectangle(MPI.COMM_WORLD, [corner_1, corner_2], [nx, ny])
    P1 = element("Lagrange", msh.basix_cell(), 1)
    ME = functionspace(msh, mixed_element([P1, P1]))
    cspace = functionspace(msh, P1)
    return msh, ME, cspace


def create_domain_urban(params, boundary_coords):
    lc = params.get("lc", 0.1)
    num_points = params.get("num_points", len(boundary_coords))

    gmsh.initialize()
    gmsh.model.add("polygon")

    # --- Add points and lines ---
    point_tags = [gmsh.model.occ.addPoint(x, y, 0, lc) for x, y in boundary_coords]

    line_tags = []
    for i in range(num_points):
        start = point_tags[i]
        end = point_tags[(i + 1) % num_points]
        line_tags.append(gmsh.model.occ.addLine(start, end))

    # --- Create surface from lines ---
    loop = gmsh.model.occ.addCurveLoop(line_tags)
    surface = gmsh.model.occ.addPlaneSurface([loop])
    gmsh.model.occ.synchronize()

    # --- Add physical groups ---
    gmsh.model.addPhysicalGroup(2, [surface], tag=1)
    gmsh.model.setPhysicalName(2, 1, "domain")

    gmsh.model.addPhysicalGroup(1, line_tags, tag=2)
    gmsh.model.setPhysicalName(1, 2, "boundary")

    # --- Mesh options ---
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)
    gmsh.model.mesh.generate(2)

    # --- Convert to dolfinx mesh ---
    msh, cell_tags, facet_tags = gmshio.model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=2)
    gmsh.finalize()

    # --- Define function spaces ---
    P1 = element("Lagrange", msh.basix_cell(), 1)
    ME = functionspace(msh, mixed_element([P1, P1]))
    cspace = functionspace(msh, P1)

    return msh, ME, cspace