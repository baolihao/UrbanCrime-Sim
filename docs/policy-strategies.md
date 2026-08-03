# Policing strategies

Police density `pi` is the canonical policy quantity.  The attenuation appearing
in the burglary equations is always computed as `exp(-pi)`; it is never stored or
exported under the name “police.”

## None

Sets `pi = 0` and adds no state variables.  Constant coefficient configurations
recover the no-police model used by the finite-element paper.

## Prescribed

Supplies `pi(x,t)` externally.  The snapshot profile currently implemented is

$$
\pi_{raw}=-\log\left[\tfrac12(1-\tanh(\kappa(A(t_s)-A_c)))\right].
$$

The profile may be scaled to a requested total budget.  Because this scaling
changes the original formula, outputs identify it as `budget_normalized_snapshot`.

## Optimal

For crime intensity `f = rho*A`, the nodal mass-lumped solution is

$$
\pi_i=\max(\log(f_i/\lambda),0),
\qquad \sum_i w_i\pi_i=M.
$$

The weights `w_i = integral(phi_i)` are assembled from the actual mesh.  This
replaces the uniform-node-area approximation in the research notebook and works
on irregular geometries.

## Delayed

Adds the dynamic fields `pi` and `H`.  The police PDE conserves the total police
mass under the natural no-flux condition.  The configured initial mean density
therefore sets the budget; no independent `M_target` is stored.

## Policy comparison

One no-police warm-up run produces an immutable snapshot at activation time.
Each policy starts from a separate copy of that snapshot.  Policies must not pass
mutable `npz` files to one another or depend on notebook execution order.

The Python ABM uses the same canonical police variable and burglary hazard
`A*exp(-pi)`.  Its delayed police transport discretizes
`div(H^2 grad(pi/H^2))` with conservative face fluxes and backward Euler, so the
configured initial police mean density fixes a conserved discrete budget.
