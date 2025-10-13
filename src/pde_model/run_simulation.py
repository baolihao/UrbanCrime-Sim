from .parameters import params
from .IMsolver import BurglaryPDEIMSolver

if __name__ == "__main__":
    solver = BurglaryPDEIMSolver(params)
    solver.run()
