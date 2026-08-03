# Canonical continuum model

UrbanCrime-Sim has one continuum formulation.  There is no legacy/paper switch.
The paper equations and historical notebooks are validation sources rather than
alternative implementations.

Let `A` denote attractiveness, `rho` criminal density, `pi` police density, `H`
delayed crime information, `ast` static attractiveness, and `q` the effective
criminal source.  The model is

\[
\partial_t A
= \nabla\!\cdot\!\left(\eta(x)\nabla(A-A^{st}(x))\right)
-A+\rho A e^{-\pi}+A^{st}(x),
\]

\[
\partial_t\rho
=\nabla\!\cdot\!\left(\nabla\rho-\frac{2\rho}{A}\nabla A\right)
-(\rho A-q(x))e^{-\pi}.
\]

Delayed policing adds

\[
\partial_t\pi
=\nabla\!\cdot\!\left(\nabla\pi-\frac{2\pi}{H}\nabla H\right),
\qquad
\partial_tH=\frac{\rho A e^{-\pi}-H}{\tau(x)}.
\]

The associated natural flux conditions are

\[
\eta\nabla(A-A^{st})\cdot n=0,
\quad
\left(\nabla\rho-2\rho\nabla A/A\right)\cdot n=0,
\quad
\left(\nabla\pi-2\pi\nabla H/H\right)\cdot n=0.
\]

### Relation to the published boundary condition

The M3AS paper writes `grad(A) . n = 0` in equation (2.17), and the police-response
paper repeats it in equation (3.11).  That condition is equivalent to the first
flux condition above only when `grad(ast) . n = 0`, including the constant-`ast`
experiments in both papers.

The agent-based update diffuses the dynamic attractiveness `B`, not the total
field `A`, and the continuum identification is `B = A - ast`.  The canonical
condition `eta*grad(A-ast) . n = 0` therefore enforces no flux of the quantity
that actually diffuses in the ABM.  For spatially varying `ast` this is a stated
heterogeneous-flux extension of the published constant-background boundary
condition, not an algebraic restatement of it.

## Backward-Euler weak form

For test functions `a`, `r`, `p`, and `h`, the implementation assembles

\[
(A-A^-,a)+\Delta t\,(\eta\nabla(A-A^{st}),\nabla a)
+\Delta t\,(A-\rho A e^{-\pi}-A^{st},a)=0,
\]

\[
(\rho-\rho^-,r)+\Delta t\,(\nabla\rho-2\rho\nabla A/A,\nabla r)
+\Delta t\,((\rho A-q)e^{-\pi},r)=0.
\]

The delayed fields use the analogous consistent weak form.  A partitioned solve
may update `H` with a documented mass-lumped pointwise backward-Euler step.

## Parameter conventions

- `q` replaces the historical names `Bbar`, `GT`, `GTS`, and `GTfactor`.
- The PDE layer does not separately expose `Gamma`, `theta`, `Sigma`, or `omega`.
- `eta > 0`, `tau > 0`, `q >= 0`, and valid states require `A > 0`, `H > 0`,
  `rho >= 0`, and `pi >= 0`.
- Scalars, analytic profiles, and loaded fields all enter the same UFL form as
  coefficients.
