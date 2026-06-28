# `src/effective/` — Effective interface parameters

Material configuration and the derived effective jump-condition parameters.

| File | Purpose |
|---|---|
| `config_parameters.m` | material + numerical parameters for the four regimes (`soft`, `intermediate`, `critical`, `stiff`); computes the effective parameters |
| `effective_params.m` | the dimensionally-correct, FE-validated effective parameters |
| `acoustic_tensor.m` | acoustic-tensor helper (also linked from `src/elements/`) |

**Effective parameters (validated against the resolved model — see `docs/FE_VALIDATION_REPORT.md`):**

```
K_n   = (λ̃ + 2μ̃) / ε_SEI      [Pa/m]   normal interface stiffness
K_t   =  μ̃        / ε_SEI      [Pa/m]   tangential interface stiffness
R_eff =  ε_SEI / D_SEI          [s/m]    diffusion resistance
β_eff =  β_SEI · ε_SEI                   chemo-mechanical coupling
G_c   =  ε_SEI · Σ_i f_i g_v,i  [J/m²]   fracture energy
δ_c   =  √(2 G_c / K_n)         [m]      critical separation
```

The normal/tangential anisotropy `K_n/K_t = 2(1−ν)/(1−2ν)` is the result the
2-D validation confirms and the 1-D model cannot test.
