# Discrete agent-based model

The Python ABM is a stochastic lattice model, not a finite-difference solver for
the continuum equations. The arrays `criminals` and `police_agents` contain
integer counts at every site. Its YAML inputs and saved fields use the
nondimensional variables of the papers.

At site $s$, a criminal commits a burglary with probability

$$
p_s=1-\exp\{-A_s e^{-\pi_s}\Delta t\}.
$$

Burglars who do not commit an event move to a neighboring site with probability
proportional to the neighboring attractiveness. New criminals are generated
with rate $\Gamma(1-\Sigma)e^{-\pi_s}$, where the configured
`generation_ratio` is $\Gamma\theta/\omega^2$. Dynamic attractiveness follows

$$
B_s^{n+1}=\left[B_s^n+\frac{\eta}{4}
\left(\sum_{r\sim s}B_r^n-4B_s^n\right)\right](1-\Delta t)
+\frac{\theta}{\omega}E_s^n.
$$

Here $E_s^n$ is the sampled integer number of burglaries. For delayed policing,
the information signal is updated explicitly from those realized events,

$$
H_s^{n+1}=H_s^n+\frac{\Delta t}{\tau}(S_s^n-H_s^n).
$$

Each police agent remains at its site with probability
$1-\exp(-\Sigma n_s p_s)$; every remaining agent moves to a neighbor with
probability proportional to the neighbor's old $H$. If all four neighboring
values are zero, the implementation uses a uniform move. This avoids the
undefined $0/0$ normalization while preserving the integer police budget
exactly.

## Discrete-to-continuum scaling

For grid spacing $h$ and nondimensional step $\Delta t$, one criminal, one
police agent, and one realized event correspond respectively to

$$
\rho_{\rm unit}=\frac{4\theta\Delta t}{\omega h^2},\qquad
\pi_{\rm unit}=\frac{4\beta\omega}{D h^2},\qquad
S_{\rm unit}=\frac{4\theta}{\omega h^2}.
$$

The initial continuous densities are converted into nonnegative integer arrays
with a fixed total equal to the ceiling of the requested population. The
fractional remainder is placed randomly without replacement; a fixed seed makes
this initialization reproducible.

`trajectory.npz` stores selected field snapshots, raw counts, the expected
continuum rate `expected_S = rho*A*exp(-pi)`, and the event-based
`realized_S`. `metrics.npz` stores lightweight spatial averages at
`time.save_every`; `summary.json` records cumulative events independently of
the output frequency.

## Implementation conventions

The maintained Python and MATLAB implementations use the following explicit
conventions:

- dimensional and nondimensional time are related by
  $t=\widetilde t/\omega$, including the conversion of `delta_t`, final time,
  and $\tau$;
- zero neighboring information produces a uniform police move instead of NaNs;
- the busy-police probability uses $n_sp_s$ as in equation (2.12);
- event totals are accumulated every step and therefore do not change with the
  snapshot frequency;
- `generation_ratio`, $\theta$, $\beta$, $D$, $\omega$, $\Sigma$, all initial
  fields, and every plotting time live in YAML rather than script constants.

The full paper configurations are under `configs/abm/`. The shorter
`delayed_square.yaml` is an inexpensive end-to-end example.

The Case 8 and Case 9 configurations use the extension paper's dimensional
value `theta: 0.2339`.
