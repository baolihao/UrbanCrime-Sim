# Reproducibility map

This file maps the **final paper figure numbers** to repository configurations
and literal commands.  Run every command from the repository root after
activating the `urbancrime` conda environment.  A row appears under "verified
entry points" only when the configuration contains the published parameter
regime and the command requests the times and fields shown in that figure.

## M3AS finite-element paper

### PDE figures: verified entry points

| Paper result | Configuration | Simulation command |
|---|---|---|
| Figure 4, Case 1 | `configs/simulations/no_police_case1_square.yaml` | `python scripts/run_pde.py configs/simulations/no_police_case1_square.yaml --snapshot 0 --snapshot 2 --snapshot 4 --snapshot 10 --snapshot 200` |
| Figure 7, Case 2 | `configs/simulations/no_police_case2_square.yaml` | `python scripts/run_pde.py configs/simulations/no_police_case2_square.yaml --snapshot 0 --snapshot 10 --snapshot 20 --snapshot 30 --snapshot 200` |
| Figure 8, Case 2 without initial noise | `configs/simulations/no_police_case2_no_noise_square.yaml` | `python scripts/run_pde.py configs/simulations/no_police_case2_no_noise_square.yaml --snapshot 0 --snapshot 110 --snapshot 120 --snapshot 130 --snapshot 200` |
| Figure 10, Case 3 | `configs/simulations/no_police_case3_square.yaml` | `python scripts/run_pde.py configs/simulations/no_police_case3_square.yaml --snapshot 0 --snapshot 2 --snapshot 4 --snapshot 10 --snapshot 200` |
| Figure 11, Case 3 without initial noise | `configs/simulations/no_police_case3_no_noise_square.yaml` | `python scripts/run_pde.py configs/simulations/no_police_case3_no_noise_square.yaml --snapshot 0 --snapshot 25 --snapshot 27 --snapshot 29 --snapshot 200` |
| Figure 14, layered `eta` | `configs/simulations/no_police_layered_eta_square.yaml` | `python scripts/run_pde.py configs/simulations/no_police_layered_eta_square.yaml --snapshot 0 --snapshot 5 --snapshot 8 --snapshot 15 --snapshot 200` |
| Figure 16, idealized square highway | `configs/simulations/no_police_highway_square.yaml` | `python scripts/run_pde.py configs/simulations/no_police_highway_square.yaml --snapshot 0 --snapshot 3 --snapshot 6 --snapshot 9 --snapshot 200` |
| Figure 17, Chicago without I-90 | `configs/simulations/no_police_chicago.yaml` | `python scripts/run_pde.py configs/simulations/no_police_chicago.yaml --snapshot 0 --snapshot 2 --snapshot 3 --snapshot 5 --snapshot 200` |
| Figure 18, Chicago with I-90 | `configs/simulations/no_police_chicago_i90.yaml` | `python scripts/run_pde.py configs/simulations/no_police_chicago_i90.yaml --snapshot 0 --snapshot 2 --snapshot 3 --snapshot 5 --snapshot 200` |

Generate the corresponding two-row field panels with:

```bash
python scripts/make_figures.py runs/no-police-case1-square --snapshots t0 t2 t4 t10 t200 --fields A rho --output runs/no-police-case1-square/figure4.png
python scripts/make_figures.py runs/no-police-case2-square --snapshots t0 t10 t20 t30 t200 --fields A rho --output runs/no-police-case2-square/figure7.png
python scripts/make_figures.py runs/no-police-case2-no-noise-square --snapshots t0 t110 t120 t130 t200 --fields A rho --output runs/no-police-case2-no-noise-square/figure8.png
python scripts/make_figures.py runs/no-police-case3-square --snapshots t0 t2 t4 t10 t200 --fields A rho --output runs/no-police-case3-square/figure10.png
python scripts/make_figures.py runs/no-police-case3-no-noise-square --snapshots t0 t25 t27 t29 t200 --fields A rho --output runs/no-police-case3-no-noise-square/figure11.png
python scripts/make_figures.py runs/no-police-layered-eta-square --snapshots t0 t5 t8 t15 t200 --fields A rho --output runs/no-police-layered-eta-square/figure14.png
python scripts/make_figures.py runs/no-police-highway-square --snapshots t0 t3 t6 t9 t200 --fields A rho --output runs/no-police-highway-square/figure16.png
python scripts/make_figures.py runs/no-police-chicago --snapshots t0 t2 t3 t5 t200 --fields A rho --output runs/no-police-chicago/figure17.png
python scripts/make_figures.py runs/no-police-chicago-i90 --snapshots t0 t2 t3 t5 t200 --fields A rho --output runs/no-police-chicago-i90/figure18.png
```

The stochastic PDE cases use repository seed `42`.  This fixes the new runs,
but does not claim pixel equality with the historical notebook figures because
the paper did not record their original random stream and software environment.
For the no-noise unstable cases, tiny floating-point perturbations select the
eventual translated pattern, so compare onset, wavelength, and amplitude rather
than hotspot coordinates.

### ABM figures: verified entry points

| Paper result | Configuration | Simulation command | Figure command |
|---|---|---|---|
| Figure 5, Case 1 | `configs/abm/no_police_case1.yaml` | `python scripts/run_abm.py configs/abm/no_police_case1.yaml` | `python scripts/make_abm_figures.py runs/no-police-abm-case-1 --times 0 3.6 10 15 200 --fields A rho --output runs/no-police-abm-case-1/figure5.png` |
| Figure 6, Case 1 late-time attractiveness | same run | same run | `python scripts/make_abm_figures.py runs/no-police-abm-case-1 --times 225 250 275 300 --fields A --output runs/no-police-abm-case-1/figure6.png` |
| Figure 9, Case 2 | `configs/abm/no_police_case2.yaml` | `python scripts/run_abm.py configs/abm/no_police_case2.yaml` | `python scripts/make_abm_figures.py runs/no-police-abm-case-2 --times 0 5 10 20 200 --fields A rho --output runs/no-police-abm-case-2/figure9.png` |
| Figure 12, Case 3 | `configs/abm/no_police_case3.yaml` | `python scripts/run_abm.py configs/abm/no_police_case3.yaml` | `python scripts/make_abm_figures.py runs/no-police-abm-case-3 --times 0 1.6 10 20 200 --fields A rho --output runs/no-police-abm-case-3/figure12.png` |

The historical MATLAB implementation remains in `legacy/matlab_abm/` as the
provenance source.  The Python commands above are the maintained executable path.

### M3AS figures not yet exposed by one command

- Figure 13 is the hotspot diameter/count sweep over `eta`; it is not the
  layered-`eta` experiment.  The individual PDE capability exists, but the sweep
  and hotspot-measurement driver still need to be packaged.
- Figure 15 reports partitioned-solver iteration counts for Case 3 without
  noise at three tolerances.  The solver runs, but iteration history is not yet
  exported in the form needed to regenerate the plot.

## Police-response extension (arXiv:2605.17709v1)

### Verified entry points

| Paper result | Configuration | Simulation command | Figure command or output |
|---|---|---|---|
| Figure 3, Hopf phase diagrams | `configs/studies/stability_figure3.yaml` | `python scripts/run_stability.py configs/studies/stability_figure3.yaml` | `runs/stability-figure-3/phase_diagrams.png` and `.npz` |
| Figure 5, Case 1 delayed-police ABM | `configs/abm/delayed_case1.yaml` | `python scripts/run_abm.py configs/abm/delayed_case1.yaml` | `python scripts/make_abm_figures.py runs/delayed-abm-case-1 --times 150 --fields A H rho pi --output runs/delayed-abm-case-1/figure5.png` |
| Figure 6, Case 2 delayed-police ABM spatial averages | `configs/abm/delayed_case2.yaml` | `python scripts/run_abm.py configs/abm/delayed_case2.yaml` | `python scripts/make_abm_figures.py runs/delayed-abm-case-2 --metrics mean_A mean_rho mean_H mean_pi --overlay --output runs/delayed-abm-case-2/figure6.png` |
| Figure 18, Case 8 delayed-police ABM | `configs/abm/delayed_case8.yaml` | `python scripts/run_abm.py configs/abm/delayed_case8.yaml` | `python scripts/make_abm_figures.py runs/delayed-abm-case-8 --times 765 768 771 774 777 780 --fields expected_S --output runs/delayed-abm-case-8/figure18.png` |
| Figure 22, Case 9 delayed-police ABM | `configs/abm/delayed_case9.yaml` | `python scripts/run_abm.py configs/abm/delayed_case9.yaml` | `python scripts/make_abm_figures.py runs/delayed-abm-case-9 --times 744 747 750 753 756 759 --fields expected_S --output runs/delayed-abm-case-9/figure22.png` |

The final extension-paper Case numbers follow Table 1 in the paper.  Several
historical notebooks in `Desktop/Crime Model/Crime-Police` use older Case names;
their parameters and saved times were checked, but those filenames are not used
as the authoritative numbering source.

### Extension studies that are implemented but not yet exact paper figures

- `configs/simulations/delayed_square.yaml` is a reusable Case-3-regime run.  It
  is not by itself an exact command for one numbered multipanel figure.
- `configs/studies/policy_comparison.yaml` contains the parameters for Figure 26
  (`eta=0.15`, activation at `t_s=200`, police mass `M=50`, fixed-policy
  `mu=5`, `A_c=1.5`).  The current plot contains post-activation averages and
  maxima, but not the no-police prefix and field insets of the published layout,
  so it is not listed above as a completed Figure 26 reproduction.
- The extension PDE field panels, spectra, return maps, and the remaining Table
  1 parameter sweeps still need dedicated figure drivers.  Generic solver runs
  are not labeled as reproductions until those analysis steps exist.

## Verification commands

```bash
pytest tests -q
python scripts/validate_config.py configs/simulations/no_police_case2_no_noise_square.yaml
python scripts/validate_config.py configs/simulations/no_police_case3_no_noise_square.yaml
python scripts/validate_config.py configs/simulations/no_police_highway_square.yaml
python scripts/validate_config.py configs/simulations/no_police_chicago_i90.yaml
```

Every completed run records the resolved YAML, Git revision, environment, mesh
summary, random seeds, and output checksums.  Inexpensive end-to-end checks live
in `configs/smoke/`; the tagged old/new numerical comparison lives in
`configs/regression/legacy_no_police_constant.yaml`.
