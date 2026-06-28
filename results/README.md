# `results/` — solver output (regenerated)

This folder holds the MATLAB `.mat` files written by the validation suite
(`case1_steady`, `case2_transient`, `case3_fracture`, `case4_convergence`,
`summary_all_cases`). They are **reproducible artifacts**, not source, and are
therefore git-ignored.

Regenerate them by running, from the repository root:

```matlab
run_all          % all phases (1-D + 2-D)
run_all('1d')    % 1-D cases 1-4 only
run_all('2d')    % 2-D anisotropy, convergence, bilayer, curvature
```
