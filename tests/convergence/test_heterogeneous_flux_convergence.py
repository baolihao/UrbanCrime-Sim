import numpy as np
import pytest

from urbancrime.models.burglary import attractiveness_flux


@pytest.mark.fenics
@pytest.mark.slow
def test_variable_eta_flux_has_second_order_l2_convergence() -> None:
    """MMS guard against replacing div(eta grad(B)) with eta*Laplacian(B)."""

    import ufl
    from dolfinx import fem, mesh
    from dolfinx.fem.petsc import LinearProblem
    from mpi4py import MPI

    errors: list[float] = []
    sizes = (8, 16, 32)
    for cells_per_side in sizes:
        domain = mesh.create_unit_square(MPI.COMM_WORLD, cells_per_side, cells_per_side)
        V = fem.functionspace(domain, ("Lagrange", 1))
        x = ufl.SpatialCoordinate(domain)
        exact = ufl.sin(ufl.pi * x[0]) * ufl.sin(ufl.pi * x[1])
        eta = 2.0 + x[0]
        forcing = (
            2.0 * ufl.pi**2 * eta * exact
            - ufl.pi * ufl.cos(ufl.pi * x[0]) * ufl.sin(ufl.pi * x[1])
        )
        trial, test = ufl.TrialFunction(V), ufl.TestFunction(V)
        bilinear = ufl.inner(attractiveness_flux(trial, 0.0, eta), ufl.grad(test)) * ufl.dx
        linear = ufl.inner(forcing, test) * ufl.dx

        boundary_facets = mesh.locate_entities_boundary(
            domain,
            domain.topology.dim - 1,
            lambda points: np.full(points.shape[1], True),
        )
        boundary_dofs = fem.locate_dofs_topological(
            V, domain.topology.dim - 1, boundary_facets
        )
        exact_boundary = fem.Function(V)
        exact_boundary.interpolate(
            lambda points: np.sin(np.pi * points[0]) * np.sin(np.pi * points[1])
        )
        boundary_condition = fem.dirichletbc(exact_boundary, boundary_dofs)
        solution = LinearProblem(
            bilinear,
            linear,
            bcs=[boundary_condition],
            petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        ).solve()
        local_error = fem.assemble_scalar(fem.form((solution - exact) ** 2 * ufl.dx))
        errors.append(float(np.sqrt(domain.comm.allreduce(local_error, op=MPI.SUM))))

    rates = [
        np.log(errors[index] / errors[index + 1])
        / np.log(sizes[index + 1] / sizes[index])
        for index in range(len(errors) - 1)
    ]
    assert min(rates) > 1.8, f"L2 errors={errors}, rates={rates}"
