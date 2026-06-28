# `src/postproc/` — Post-processing and figures

Error metrics and figure generation for the validation.

| File | Purpose |
|---|---|
| `compute_l2_error.m` | relative `L²` error between two 1-D fields (trapezoidal norm; handles the doubled interface node) |
| `plot_field_comparison.m` | FULL vs EFF field comparison (displacement / concentration / stress) |
| `plot_cohesive_curves.m` | traction–separation curves at several average concentrations |
| `plot_convergence_study.m` | log–log convergence plot (Figure 5: rate vs `ε/L` for the four regimes) |

Figures are written to `programmes/figures/`; `.mat` results to `programmes/results/`.
