# `tests/` — Validation cases

Each test is runnable standalone (`matlab -batch test_name`) — they add the
source tree to the path automatically — or together via `run_all.m`.

## 1-D validation (Section 6.2–6.5 of the paper)

| File | Paper § | Checks | Expected |
|---|---|---|---|
| `main_validation.m` | 6 | driver: runs cases 1–4, writes figures + `results/` | global PASS |
| `test_case1_steady.m` | 6.2 | steady trilayer; FULL vs EFF | u **0.00 %**, c **1.55 %** |
| `test_case2_transient.m` | 6.3 | transient cyclic loading | **0.01 %** |
| `test_case3_fracture.m` | 6.4 | cohesive law; `Gc` three independent ways | agree to **4 %** |
| `test_case4_convergence.m` | 6.5 | `O(ε)` convergence, four regimes | rates **1.67–1.91** |
| `test_scaled_modulus.m` | 6.5 | rates independent of the modulus scaling `C^ε=(ε/L)^α C̃` | identical curves |

## 2-D validation (plane strain — the part 1-D cannot test)

| File | Checks | Expected |
|---|---|---|
| `test_2d_anisotropy.m` | normal/tangential interface stiffness | `K_n=(λ+2μ)/ε`, `K_t=μ/ε`, `K_n/K_t = 2(1−ν)/(1−2ν) = 3.27`; an isotropic interface fails in shear (jump ×0.31) |
| `test_2d_convergence.m` | normal **and** shear interface-jump convergence | rates **0.97 / 0.97** ≈ `O(ε)` |
| `test_2d_bilayer.m` | bilayer (inorganic+organic) series interface law | `K_n,K_t` match series formula exactly; `K_n` soft-dominated (0.18× the naive average); `R_eff` blocking-dominated |
| `test_2d_curved.m` | curved interface (axisymmetric radial) | error ∝ `ε/R` (slope 1.00); flat limit `R→∞` recovers the flat interface (err 3e-4); ≤0.75% for real particles `R≥20 µm` |
| `make_profiles_2d.m` | solution-profile snapshots (resolved vs homogenized) | writes `figures/fig_profiles_2d.png`, `fig_fields_2d_*.png` |

All results come from actual MATLAB finite-element runs (no commercial FE).
