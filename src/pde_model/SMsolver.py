import numpy as np
from dolfinx import fem, log
from dolfinx.fem import Function
from numpy import linalg as LA
from dolfinx.fem.petsc import LinearProblem
from dolfinx.mesh import CellType
from ufl import dx, grad, inner, div
from petsc4py import PETSc
import ufl

from .mesh_utils import create_domain

class BurglaryPDESMSolver:
    def __init__(self, params,Anoise, rhonoise, msh, ME, cspace, A0_fun=None, Bbar_fun=None):
        self.params = params
        self.msh, self.ME, self.cspace = msh, ME, cspace
        self.dt = params["dt"]
        self.T = params["T"]
        self.Nx, self.Ny = params["Nx"],params["Ny"]
        self.A_st = params["A_st"]
        self.rho0 = params["rho0"]
        self.max_Iter = params["max_Iter"]
        # Load initial conditions
        Ai = 0.1* np.random.rand((self.Nx+1)* (self.Ny+1)) + self.A_st
        rhoi = 0.01* np.random.rand((self.Nx+1)* (self.Ny+1)) + self.rho0

        # Initialize functions
        U0, dofsU = self.ME.sub(0).collapse()
        V0, dofsV = self.ME.sub(1).collapse()
        self.A0_fun = A0_fun
        self.Bbar_fun = Bbar_fun
        self.u0 = Function(U0)
        self.v0 = Function(V0)# solution from previous converged step
        self.uB = Function(U0)
        self.vrho = Function(V0)# solution from previous converged step
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

    def _initialize_fields(self, Ai, rhoi):
        self.u0.x.array[:] = Ai + self.Bbar.x.array[:] 
        self.v0.x.array[:] = rhoi

        self.uB.x.array[:] = Ai + self.Bbar.x.array[:] 
        self.vrho.x.array[:] = rhoi
    
        self.u0.x.scatter_forward()
        self.v0.x.scatter_forward()
        self.A0.x.scatter_forward()
        self.uB.x.scatter_forward()
        self.vrho.x.scatter_forward()


    def _define_problem(self):
        p = ufl.TestFunction(self.cspace)
        q = ufl.TestFunction(self.cspace)
        u = ufl.TrialFunction(self.cspace)
        v = ufl.TrialFunction(self.cspace)# current solution
      
        dt, eta, Bbar = self.dt, self.params["eta"], self.params["Bbar"]
        F0_left = inner(u,p) * dx  + eta * dt * inner(grad(u), grad(p)) * dx  -   dt * inner(u*self.vrho-u, p) * dx
        F0_right = inner(self.u0, p) * dx - dt * inner(eta*div(grad(self.A0))- self.A0, p) * dx
        F1_left = (inner(v, q) * dx
                    + dt * inner(grad(v), grad(q)) * dx
                    - 2 * dt * inner(v / self.uB * grad(self.uB), grad(q)) * dx
                    + dt * inner(self.uB * v, q) * dx)
        F1_right = inner(self.v0, q) * dx + dt * inner(Bbar, q) * dx 
        self.F0_left, self.F0_right = F0_left, F0_right
        self.F1_left, self.F1_right = F1_left, F1_right


    def run(self):
        t, step = 0.0, 0
        while t < self.T:
            t += self.dt
            k = 1
            step += 1
            while k < self.max_Iter+1:
                problem1 = LinearProblem(self.F0_left, self.F0_right, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
                u_new = problem1.solve()
                Bnorm = LA.norm(np.subtract(u_new.x.array[:],self.uB.x.array[:]),2)/LA.norm(self.uB.x.array[:],2)
        
                self.uB.x.array[:] = u_new.x.array[:]      
                problem2 = LinearProblem(self.F1_left, self.F1_right, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
                v_new = problem2.solve()


                rhonorm = LA.norm(np.subtract(v_new.x.array[:],self.vrho.x.array[:]),2)/LA.norm(self.vrho.x.array[:],2)
        

                self.vrho.x.array[:] = v_new.x.array[:]

        
        
                if (Bnorm <1e-9 and  rhonorm <1e-9):
                    self.u0.x.array[:] = u_new.x.array[:]
                    self.v0.x.array[:] = v_new.x.array[:]
                    break
            
                k+=1

            print(f"Step {step}: num iterations: {k}")
        print("Simulation finished.")
