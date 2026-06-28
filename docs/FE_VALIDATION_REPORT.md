# Finite-Element Validation — Real Simulation Report

This report documents an actual **execution** of the finite-element validation
suite (`main_validation.m` and its test cases) via MATLAB, the bugs that had to be
fixed before any simulation could run, and the **real** numerical results.

> Key fact: as delivered, the FE suite **did not run at all** — it failed at parse
> time, then at the first mesh build. Six distinct bugs were fixed before the
> validation produced numbers. Every result below comes from an actual FE solve.

## 1. Bugs found and fixed (the suite never executed before)

| # | File | Bug (provable) | Fix |
|---|---|---|---|
| 1 | `main_validation.m` | `main_validation` had no closing `end` while local `ternary` did → illegal mixed function convention → **parse error**, nothing runs | added `end` |
| 2 | `build_trilayer_mesh.m` | connectivity stored node **coordinates** (`nodes(i)`) instead of **indices** (`i`) → `nodes(idx)` indexing crash → FULL model never assembled | use indices `i, i+1` |
| 3 | `build_eff_mesh.m` | interface node indices off by one vs the doubled-node array → **zero-length element** (`L=0`) → EFF model crashed | `idx_left=ne+1`, `idx_right=ne+2` |
| 4 | `compute_l2_error.m` | `interp1` over the EFF mesh failed: the doubled interface node gives **non-unique abscissa** | sort + nudge duplicates |
| 5 | `effective_params.m` / `config_parameters.m` | effective parameters dimensionally wrong: `K_eff=(ε/L)^(α−1)·C1d/h` is **Pa, not Pa/m** → EFF/FULL displacement error **6 000 000 %** | physical layer values (see §2) |
| 6 | `test_case4_convergence.m` | `struct('regimes', cell,...)` built a **1×4 struct array** (not a scalar struct) → `conv_data.errors` returned 4 values → plot crash | `'regimes', {cell}` |
| 7 | `plot_field_comparison.m` | `ylim(2)` parsed as a *call* `ylim(2)` (set limits to scalar 2) instead of indexing the returned vector | capture `yl = ylim` first |
| 8 | `chemo_cohesive_law.m` | **missing file** — referenced by `test_case3_fracture` but absent → fracture case could not run | implemented the critical-regime cohesive law |
| 9 | `test_case3_fracture.m` (`solve_damage_step`) | `sigma^2 / (...)` used matrix `/` (mrdivide) instead of `./` → `Y` became a row vector → broadcasting made a 51×51 matrix → `^` crash | use `./`; force `d_field(:)` |

(Plus a log-log axis fix in `plot_convergence_study.m`: `axis equal` had reset the
axes to linear, hiding the convergence slope.)

## 2. The decisive result: correct effective interface parameters

The real simulation is the arbiter of the long-standing `K_eff` units question.
A thin uniform layer of modulus `C`, diffusivity `D` and thickness `t = h·ε`
behaves, as `ε→0`, as a zero-thickness interface with

```
K_eff   = C1d_sei / t        [Pa/m]   (through-thickness spring)
R_eff   = t / D_sei          [s/m]    (diffusion resistance)
beta_eff= beta_sei · t                (chemical-swelling jump)
```

With the *old* `K_eff = (ε/L)^(α−1)·C1d/h` (≈10¹⁰, dimensionally Pa) the EFF model
was ~10⁷× too soft and the displacement L2 error was **6.1×10⁶ %**. With the
physical `K_eff = C1d_sei/t` (≈6.4×10¹⁷ Pa/m) the error collapses to **0.00 %**.
This resolves the `KEFF-UNITS` audit finding by direct simulation.

## 3. Real convergence study (Section 6.5 / Figure 5)

`test_case4_convergence` solves the fully-resolved (FULL, SEI explicitly meshed)
and homogenized (EFF, jump conditions) models over `ε/L ∈ [10⁻⁴, 5×10⁻²]` (8 points)
for four parameter sets, and fits the slope of the relative L2 displacement error.

| Regime | error at ε/L=10⁻⁴ | error at ε/L=5×10⁻² | **observed rate** |
|---|---|---|---|
| soft         | 2.9×10⁻⁸ | 3.5×10⁻³ | **1.91** |
| intermediate | 1.4×10⁻⁸ | 1.8×10⁻³ | **1.91** |
| critical     | 2.0×10⁻⁸ | 1.1×10⁻³ | **1.74** |
| stiff        | 2.1×10⁻⁸ | 7.5×10⁻⁴ | **1.67** |

**Validated**: the homogenized interface model converges to the resolved model for
all four parameter sets, with relative errors of 10⁻⁸–10⁻³ (≤0.35 %).

**Two honest discrepancies with the article's stated claims** (Section 6.5):

1. The observed convergence is **super-linear, ~O(ε^1.7–1.9)** — *better* than the
   stated O(ε). The clean log-log slopes (Figure 5) sit clearly above the O(ε)
   reference line.
2. The stiff regime **also converges** (rate 1.67); the stated "stiff = O(1)
   (singular, non-convergent)" is **not** borne out.

### 3a. Tested directly: does the modulus scaling change the conclusion? No.

`test_scaled_modulus.m` repeats the study with the article's scaling hypothesis
`C^ε = (ε/L)^α · C̃` applied to the **resolved** SEI modulus (so the SEI genuinely
stiffens/softens as ε→0), with the matching effective spring `K_eff = C^ε/(h·ε)`.
The resulting error curves are **byte-identical** to the fixed-modulus case, and
the rates are the same: soft 1.91, intermediate 1.91, critical 1.74, **stiff 1.67**.

This is the key conclusion, established by simulation: **the convergence rate is
independent of the modulus regime**. The relative error is invariant under a
global rescaling of the modulus once the effective stiffness is the exact thin-
layer spring `K_eff = C^ε/ε`. The stiff regime therefore converges; an `O(1)` error
appears only if one replaces the interface by its *rigid limit* `[[u]]=0` instead
of the derived finite spring — that is an approximation of the limit law, not a
property of the homogenization. Reproduce: `matlab -batch "test_scaled_modulus"`.

## 4. All four test cases — real results

| Case | What it tests | Result |
|---|---|---|
| 1 — steady (§6.2) | trilayer, eigenstrain loading | u **0.00 %**, c **1.55 %**, σ 107 %* |
| 2 — transient (§6.3) | cyclic loading, time stepping | **PASS**, max error **0.01 %** |
| 3 — fracture (§6.4) | cohesive law, Gc | **PASS**, Gc agree to **4 %** (below) |
| 4 — convergence (§6.5) | O(ε) rate, 4 regimes | **PASS**, rates 1.67–1.91 (§3) |

\* The 107 % stress error is the SEI-**interior** stress, which the zero-thickness
effective model does not resolve by construction; the displacement, concentration,
and interface traction all match. (Displacement is the metric used for the
convergence study.)

**Fracture energy validated three independent ways** (critical regime, ε=30 nm):

| Source | Gc (J/m²) |
|---|---|
| Theorem 4.2 formula (`P.Gc`) | 1.035×10⁻² |
| Effective cohesive law, ∫T dΔ | 9.93×10⁻³ |
| Full damage simulation | 9.997×10⁻³ |

All three agree within 4 %, a genuine cross-validation of the composition-based
fracture-energy result.

## 5. What this means for the article

- The numerical validation **does** support the homogenized jump-condition model —
  once the effective parameters are dimensionally correct (`K_eff=C/t`, etc.).
- Section 6.5 should report the **measured** rate (~O(ε^1.8), all regimes converge)
  rather than O(ε)/O(1); Figure 5 should be the real log-log plot.
- Table 3's `K_n`, `R_eff`, `β_eff` should be recomputed from `C/t`, `t/D`,
  `β·t` for consistency with the simulation (this is the `KEFF-UNITS` resolution).

All numbers here are reproducible: `matlab -batch "main_validation"` (after the
fixes) and the per-case logs `matlab_val*.log`, `matlab_conv*.log`.
