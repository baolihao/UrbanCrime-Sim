# UrbanCrime-Sim

UrbanCrime-Sim is a configuration-driven research codebase for continuum and
agent-based models of residential burglary and police response.

The canonical continuum formulation is the heterogeneous conservative-flux
model.  Scalar parameters are represented as constant coefficients, so the same
weak form also supports spatially varying `eta`, `ast`, `q`, and (for delayed
policing) `tau` without special scalar branches.

## Model family

The persistent fields are attractiveness `A` and criminal density `rho`.  A
policing strategy supplies police density `pi`:

- `none`: no police, `pi = 0`;
- `prescribed`: a fixed or externally supplied police field;
- `optimal`: a mass-constrained instantaneous allocation;
- `delayed`: dynamic police density `pi` coupled to delayed crime information `H`.

See [docs/model.md](docs/model.md) and
[docs/policy-strategies.md](docs/policy-strategies.md) for equations and numerical
conventions.
The discrete stochastic implementation is documented separately in
[docs/agent-based-model.md](docs/agent-based-model.md).
The no-hidden-parameter rules are documented in
[docs/configuration.md](docs/configuration.md).

## Environment

FEniCSx and its MPI/PETSc dependencies are provided through conda-forge:

```bash
conda env create -f environment.yml
conda activate urbancrime
```

`pyproject.toml` contains package metadata and the pure-Python dependencies.  It
does not attempt to install the native FEniCSx stack through pip.

## Configuration


```text
configs/
├── simulations/   # one continuum PDE run
├── abm/           # one Python agent-based run
├── studies/       # sweeps and multi-strategy comparisons
├── regression/    # small tagged-old/new numerical reference
└── smoke/         # inexpensive CI checks
```

Validate a simulation config with:

```bash
python scripts/validate_config.py configs/simulations/no_police_square.yaml
```

Run the homogeneous stability study with:

```bash
python scripts/run_stability.py configs/studies/stability_figure3.yaml
```

Run the continuum model, agent-based model, and policy comparison with:

```bash
python scripts/run_pde.py configs/simulations/no_police_square.yaml
python scripts/run_abm.py configs/abm/delayed_square.yaml
python scripts/compare_policies.py configs/studies/policy_comparison.yaml
```

The equivalent installed commands are `urbancrime-pde`, `urbancrime-abm`, and
`urbancrime-compare-policies`. Each run writes its resolved configuration, Git
revision, environment metadata, numerical output, and SHA-256 checksums beneath
`runs/`.

## Verification

The complete local suite includes pure unit tests, a tagged-old/new regression,
real DOLFINx smoke tests, time and heterogeneous-flux spatial convergence,
delayed-police mass conservation, and the stochastic ABM-to-PDE reaction-rate
limit:

```bash
pytest tests -q
```

The three notebooks in `notebooks/` are intentionally thin executable front
ends.  They contain no copied weak forms, solver implementations, or scientific
parameter dictionaries.

All generated data belong under `runs/`, which is not tracked by Git.

## Papers

The model family builds on:

- B. Hao, K. Mily, A. Quaini, and M. Zhong, “A finite element framework for
  simulating residential burglary in realistic urban geometries,” *Mathematical
  Models and Methods in Applied Sciences* 36(5), 2026.
  [DOI](https://doi.org/10.1142/S0218202526500193)
- B. Hao, K. Mily, A. Quaini, and M. Zhong, “Crime hotspot dynamics in
  residential burglary models with police response,” 2026.
  [arXiv:2605.17709](https://arxiv.org/abs/2605.17709)

Figure-to-command mappings are maintained in [docs/reproduce.md](docs/reproduce.md).
The executable MATLAB reference ABM and its smoke tests are documented in
[matlab/README.md](matlab/README.md); pre-refactor sources remain under
`legacy/matlab_abm/` for provenance only.

## License and data terms

The UrbanCrime-Sim software is released under the
[BSD 3-Clause License](LICENSE).  This license covers the software, not bundled
third-party geographic data.  The derived Chicago boundary remains subject to
the City of Chicago source terms and attribution requirements; see
[data/boundaries/NOTICE.md](data/boundaries/NOTICE.md) and
[docs/licensing.md](docs/licensing.md).
