# CrimeFEM: A Finite Element and Agent-Based Modeling Framework for Urban Crime Dynamics
<div align="center">
  <img src="demo.gif" alt="Demo" width="600">
</div>

This repository provides the code accompanying the paper as a base for multiscale, multi-physics extensions, integrating agent-based models and PDE solvers:

**A finite element framework for simulating residential burglary in realistic urban geometries**  
_Baoli Hao, Kamrun Mily, Annalisa Quaini, and Ming Zhong (2025)_  
arXiv: [2508.11055](https://arxiv.org/abs/2508.11055)


It contains:
- Two **finite element solvers (Python)** for the PDE mean-field model.
- An **agent-based model (MATLAB)** representing individual criminal motion.
- Tools for data analysis, and visualization.

---

## 0. Model Overview
We simulate the evolution of the *attractiveness field* and *criminal density* over realistic city geometries,
derived as a mean-field limit from an underlying agent-based process.
The PDE model resembles a nonlinear Keller–Segel system with Neumann-type boundary conditions.

---

## 1. Repository Structure
- `src/pde_model/`: Python FEM implementation.
- `src/agent_based_model/`: Agent-based MATLAB code.
- `demo/`: Reproducibility and visualization notebooks.
- `data/`: Input data (urban maps, parameters, etc.).

---

## 3. Installation and Usage
### 1️⃣ Requirements
Before installing, ensure that you have the following dependencies available:
**Core Requirements (Python side):**
- Python ≥ 3.10  
- [FEniCSx / dolfinx](https://fenicsproject.org) (latest stable)
- `petsc4py`
- `mpi4py`
- `gmsh` (with Python API)
- `pyvista`, `pyvistaqt` (optional, for visualization)
- `numpy`

**Optional (for MATLAB ABM):**
- MATLAB R2023a or newer
- Statistics and Parallel Computing Toolboxes (recommended)

You can install most Python dependencies with:

```bash
pip install --upgrade pip
pip install -r requirements.txt
```
### 2️⃣ Environment Setup
If you use *conda* (recommended):
```bash
conda create -n burglary python=3.10
conda activate burglary
```
Then install FEniCSx and related libraries:
```bash
pip install fenics-dolfinx mpi4py petsc4py gmsh pyvista pyvistaqt numpy
```
Alternatively, you may use *Docker or micromamba* to install FEniCSx for better compatibility.

### 3️⃣ Clone the Repository
Choose a folder where you want to install the package and clone the repository:
```bash
cd ~
git clone https://github.com/baolihao/UrbanCrime-Sim.git
cd UrbanCrime-Sim
```

### 4️⃣ Run the PDE Simulation (IPython / Jupyter)
All main Python code for the PDE solver is located under:
```bash
src/pde_model/
```
You can run the PDE solver interactively using Jupyter Notebook:
```bash
jupyter notebook
```
Then open:
```bash
demo/pde/PDE_SM_Demo.ipynb
```
and execute the notebook cells step by step.

### 5️⃣ Run the Agent-Based Model (MATLAB)
The agent-based model scripts are located in:
```bash
src/agent_based_model/
```
To run them in MATLAB:
```bash
cd src/demo/ode
run test_system.m
```

## 4. Citation

If you use this code, please cite:
```
@article{hao2025finite,
  title={A finite element framework for simulating residential burglary in realistic urban geometries},
  author={Hao, Baoli and Mily, Kamrun and Quaini, Annalisa and Zhong, Ming},
  journal={arXiv preprint arXiv:2508.11055},
  year={2025}
}
```
