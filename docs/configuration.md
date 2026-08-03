# Configuration contract

Configuration files are named by the computation they perform.  Paper names,
journal acronyms, and formulation switches do not belong in configuration paths.

Every quantity that changes a numerical result must either be present in the YAML
input or be a derived quantity recorded in `config.resolved.yaml`.  Scientific
defaults must not be hidden in notebooks, scripts, solver constructors, or plotting
code.

## Parameter groups

- Model coefficients: `eta`, `ast`, `q`, and delayed `tau`.
- Policing: strategy, activation time, budget, profile parameters, positivity
  floor, optimal-control update cadence, and root-solver controls.
- Domain and mesh: bounds or boundary file, optional polygon scale factor, cell
  type, degree, resolution, and boundary size-field parameters.
- Time: start, final time, step, and output cadence.
- Initial state: all field profiles, noise distributions, amplitudes, sparsity,
  and seeds.
- Solver: method, algebraic tolerances, state-positivity tolerance, iteration
  limit, block order, and PETSc options.
- Stability studies: homogeneous parameters, domain length, spectral cutoff,
  root bracket, sampling density, and transversality step.
- Outputs: requested fields, formats, snapshot times, and output root.

The historical identifiers `Bbar`, `GT`, `GTS`, and `GTfactor` are represented by
the single PDE coefficient `q`.  Separate `Gamma`, `theta`, `Sigma`, and `omega`
belong only to the ABM or nondimensionalization layer.

## Coefficients

YAML coefficients use one of three forms:

```yaml
eta: {kind: constant, value: 0.15}
```

```yaml
eta:
  kind: analytic
  profile: sinusoidal_x
  parameters:
    offset: 0.15
    amplitude: 0.05
    wavelength: 5.0
    phase: 0.0
```

```yaml
eta:
  kind: file
  path: data/coefficients/eta_points.npz
  field: values
```

`.npy` stores one value per local degree of freedom and is intended for serial
checkpoints on an identical mesh.  Portable `.npz` fields contain `points=(n,2)`
and a named value array; they are interpolated to the active mesh, with nearest
neighbor continuation outside the point-cloud convex hull.

Python callers may pass a callable directly.  YAML must use registered profile
names rather than executable expression strings.

## Initial conditions and noise

`ast_plus_q_plus_noise` constructs the paper convention
`A(0) = ast + q + noise`.  Noise is added to this baseline; sparsity never
removes the baseline itself.  Use `A_sparsity` and `rho_sparsity` when the two
fields use different Bernoulli masks, as in the highway experiments.  A common
`sparsity` remains available when the masks have the same density.

Polygon boundaries are stored once in normalized coordinates.  A run may set
`domain.scale_factor` to reproduce a paper coordinate system without committing
duplicate boundary arrays.
