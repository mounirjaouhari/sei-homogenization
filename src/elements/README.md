# `src/elements/` — Finite-element routines

Element-level stiffness/force routines for the coupled chemo-mechanical problem.

| File | Dim | Purpose |
|---|---|---|
| `elastic_element_stiff.m` | 1D | 2-node bar element: `Ke = (C1d/L)[1 -1; -1 1]` + chemical-eigenstrain force `C1d·β(c̄−c0)[−1;1]` |
| `diffusion_element_stiff.m` | 1D | 2-node diffusion element: conductivity `K_c` and mass `M_c` |
| `interface_element.m` | 1D | zero-thickness interface: `K_int = K_eff[1 -1;-1 1]`, swelling force `K_eff·β_eff(c̄−c0)[−1;1]` |
| `elastic_q4_2d.m` | 2D | **Q4 bilinear plane-strain** element (2×2 Gauss): 8×8 stiffness + isotropic eigenstrain force. Plane-strain `D` matrix from `(E,ν)` |
| `acoustic_tensor.m` | — | acoustic (Christoffel) tensor `Q_n = (λ+μ) n⊗n + μ I` used for the through-thickness compliance |

The 2D interface (normal `K_n` + tangential `K_t` springs) is assembled directly
in `src/solvers/solve_mech_2d.m` as node-pair springs weighted by the tributary
length along the interface.
