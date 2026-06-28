# `src/solvers/` — Assembly and solvers

Global assembly, boundary conditions and the steady / transient / 2-D solvers.

| File | Dim | Purpose |
|---|---|---|
| `assemble_full_model.m` | 1D | global `K_u, f_u, K_c, M_c` for the resolved trilayer (SEI meshed) |
| `assemble_eff_model.m` | 1D | global matrices for the effective model (bulk + zero-thickness interface) |
| `solve_steady_state.m` | 1D | steady chemo-mechanics: solve diffusion `K_c c = f_c`, then mechanics `K_u u = f_u(c)` |
| `solve_transient_cyclic.m` | 1D | transient cyclic loading (implicit Euler in time) |
| `solve_mech_2d.m` | **2D** | **plane-strain** elasticity (FULL/EFF) with anisotropic interface (`K_n`, `K_t`), optional eigenstrain and **periodic-x** option |
| `chemo_cohesive_law.m` | — | critical-regime traction–separation law `T=(1−d)²K_eff Δ_eff`, calibrated to `Gc` and `δ_c` |
| `compute_damage_field.m` | — | damage field wrapper around the cohesive law |
| `apply_dirichlet_bc.m` | — | Dirichlet boundary-condition application (penalty/elimination) |
| `safe_spmatrix.m` | — | robust sparse-matrix construction helper |

The 2D solver supports `opts.periodic_x = true` (master–slave tying of the
left/right boundaries) to realise clean uniaxial-strain and simple-shear states
for the anisotropy and convergence tests.
