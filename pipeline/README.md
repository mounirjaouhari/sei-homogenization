# `pipeline/` — Materials-Project integration pipeline (Python)

Offline-capable Python pipeline that assembles the SEI constituent data,
homogenises it to effective interphase properties, and regenerates the LaTeX
tables, the MATLAB parameter file and the figures. It reuses the verified
physics in `../../Skills/sei_physics.py` (single source of truth), located by an
upward search so the pipeline is location-independent.

| File | Purpose |
|---|---|
| `main.py` | full pipeline: fetch → homogenize → effective params → MATLAB/LaTeX/figures |
| `recompute_effective.py` | self-contained core: regenerates `data/effective_parameters_corrected.csv` and prints a before/after comparison |
| `growth_fade.py` | **calibration-free capacity fade** from a growing SEI (two-scale time): parabolic growth ε∝√N → universal √N fade + composition-dependent cracking threshold; writes `figures/fig_capacity_fade.png` |
| `sensitivity_design.py` | **robustness + design**: log-elasticities (Gc∝σ_r², cracking↔mechanics vs fade↔growth decoupled) + dimensionless design master-curve Δε_crit(f_org) + predicted-vs-measured table; writes `figures/fig_sensitivity.png`, `figures/fig_master_curve.png` |
| `example_usage.py` | four worked examples (single constituent, ASEI search, Gc heatmap, scheme comparison) |
| `src/_engine.py` | bridge to `Skills/sei_physics.py` |
| `src/config.py` | compositions, exponents, output directories, `KAPPA` |
| `src/fetch_mp_data.py` | constituent properties (offline literature; optional Materials Project) |
| `src/homogenize.py` | Mori–Tanaka (elastic), Bruggeman (diffusion), Vegard, fracture mixing |
| `src/compute_effective.py` | effective parameters per composition + sweep |
| `src/generate_matlab.py`, `src/generate_latex.py` | emit `config_parameters_mp.m`, Table 1/3, §5 numbers |
| `src/plot_composition.py` | Gc-vs-composition, phase diagram, parameter comparison |

## Run

```bash
python pipeline/recompute_effective.py     # corrected CSV + before/after
python pipeline/main.py --offline          # full pipeline (no network)
python pipeline/main.py --offline --sweep  # + Gc composition sweep (Figure 4)
python pipeline/example_usage.py           # the four examples
```

Requirements: see `requirements.txt` (numpy, pandas, scipy, matplotlib;
`mp-api` optional for online Materials-Project queries).
