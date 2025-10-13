import numpy as np
from dolfinx import fem, log
from dolfinx.fem import Function
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.mesh import CellType
from ufl import dx, grad, inner, div
from petsc4py import PETSc
import ufl

from .mesh_utils import create_domain

class BurglaryPDEIMSolver:
    def __init__(self, params, Anoise, rhonoise,msh, ME, cspace, A0_fun=None, Bbar_fun=None):
        self.params = params
        self.msh, self.ME, self.cspace = msh, ME, cspace
        self.dt = params["dt"]
        self.T = params["T"]
        self.Nx, self.Ny = params["Nx"],params["Ny"]
        self.Lx, self.Ly = params["Lx"],params["Ly"]
        self.A_st = params["A_st"]
        self.rho0 = params["rho0"]
        # Load initial conditions


        # Initialize functions
        self.A0_fun = A0_fun
        self.Bbar_fun = Bbar_fun
        self.c = Function(self.ME)
        self.c0 = Function(self.ME)
        self.A0 = Function(self.cspace)
        self.Bbar = Function(self.cspace)
        if self.A0_fun is not None:
            self.A0.interpolate(lambda x: self.A0_fun(x, self.params))
        else:
            self.A0.x.array[:] = self.params["A_st"]

        if self.Bbar_fun is not None:
            self.Bbar.interpolate(lambda x: self.Bbar_fun(x, self.params))
        else:
            self.Bbar.x.array[:] = self.params["Bbar"]

       
        Ai = Anoise + self.A0.x.array[:] 
        rhoi = rhonoise + self.rho0
        self._initialize_fields(Ai, rhoi)
        self._define_problem()
        self._configure_solver()

    def _initialize_fields(self, Ai, rhoi):
        U0, dofsU = self.ME.sub(0).collapse()
        V0, dofsV = self.ME.sub(1).collapse()
        self.c.x.array[dofsU] = Ai+self.Bbar.x.array[:]
        self.c.x.array[dofsV] = rhoi
        self.c.x.scatter_forward()
        self.A0.x.scatter_forward()
        self.Bbar.x.scatter_forward()
        self.c0.x.array[:] = self.c.x.array


    def _define_problem(self):
        p, q = ufl.TestFunctions(self.ME)
        u, v = ufl.split(self.c)
        u0, v0 = ufl.split(self.c0)
        dt, eta = self.dt, self.params["eta"]

        F0 = (inner(u, p) - inner(u0, p)
              + eta * dt * inner(grad(u), grad(p))
              + dt * inner(u, p)
              - dt * inner(u * v, p)
              + dt * inner(eta * div(grad(self.A0)) - self.A0, p)) * dx

        F1 = (inner(v, q) - inner(v0, q)
              + dt * inner(grad(v), grad(q))
              - 2 * dt * inner(v / u * grad(u), grad(q))
              + dt * inner(u * v - self.Bbar, q)) * dx

        self.problem = NonlinearProblem(F0 + F1, self.c)

    def _configure_solver(self):
        self.solver = NewtonSolver(self.msh.comm, self.problem)
        self.solver.convergence_criterion = "incremental"
        self.solver.rtol = np.sqrt(np.finfo(float).eps) * 1e-6
        
        ksp = self.solver.krylov_solver
        opts = PETSc.Options()
        prefix = ksp.getOptionsPrefix()
        opts[f"{prefix}ksp_type"] = "preonly"
        opts[f"{prefix}pc_type"] = "lu"
        ksp.setFromOptions()

    def run(self):
        t, step = 0.0, 0
        while t < self.T:
            t += self.dt
            step += 1
            r = self.solver.solve(self.c)
            print(f"Step {step}: num iterations: {r[0]}")
            self.c0.x.array[:] = self.c.x.array
        print("Simulation finished.")
