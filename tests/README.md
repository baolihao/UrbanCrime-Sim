# Tests

- `unit/`: pure configuration, discrete-ABM, policy, stability, and utility tests.
- `regression/`: tagged-old/new agreement plus real two-field and four-field DOLFINx smoke runs.
- `convergence/`: backward-Euler time convergence and a variable-`eta` manufactured
  solution with second-order P1 `L2` convergence.
- `validation/`: stochastic ABM burglary rates against the PDE reaction limit.

Tests requiring FEniCSx use the `fenics` marker and are skipped when the native
conda environment is unavailable; long tests use `slow`.

The maintained MATLAB ABM has a separate fast suite:

```matlab
addpath(genpath('matlab'));
results = runtests('matlab/tests');
assertSuccess(results);
```
