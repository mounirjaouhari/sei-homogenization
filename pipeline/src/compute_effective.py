"""
compute_effective.py — effective interface parameters per SEI composition.

All power laws use the dimensionless base eta = eps/L (corrected convention);
Gc uses the volumetric g_v (Option B). Physics delegated to Skills/sei_physics.
"""

from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd

from _engine import P
from config import (SEI_COMPOSITIONS, SEI_THICKNESSES, REGIME_EXPONENTS,
                    L_MACRO, E_ELECTRODE)
from homogenize import (homogenize_sei_elastic, homogenize_sei_diffusivity,
                        homogenize_sei_beta, homogenize_fracture_energy)


def identify_regime(E_SEI: float, eps: float, gamma: float = 1.0) -> Dict:
    """
    Identify alpha from the homogenized modulus with the dimensionless base
    eta = eps/L (corrected formula), and classify the mechanical/diffusion regime.
    """
    alpha = P.identify_alpha(E_SEI, eps, L=L_MACRO, E_tilde=E_ELECTRODE)
    return {
        "alpha": alpha,
        "mechanical_regime": P.classify_mechanical_regime(alpha)["regime"],
        "diffusion_regime": P.classify_diffusion_regime(gamma)["regime"],
        "well_posed": P.well_posed(alpha, gamma),
    }


def compute_effective_parameters(comp_name: str, df: pd.DataFrame) -> Dict:
    """Full effective-parameter set for one named SEI composition."""
    comp = SEI_COMPOSITIONS[comp_name]
    eps = SEI_THICKNESSES[comp_name]
    exps = REGIME_EXPONENTS[comp_name]
    gamma, delta = exps["gamma"], exps["delta"]

    el = homogenize_sei_elastic(df, comp)
    D_sei = homogenize_sei_diffusivity(df, comp)
    beta_sei = homogenize_sei_beta(df, comp)

    reg = identify_regime(el["E"], eps, gamma=gamma)
    alpha = reg["alpha"]

    eff = P.effective_parameters(el["E"], el["nu"], D_sei, beta_sei,
                                 eps, alpha, gamma, delta, L=L_MACRO, h=1.0)

    g_v_sei = homogenize_fracture_energy(df, comp)   # J/m^3
    Gc = eps * g_v_sei                               # J/m^2
    delta_c = P.critical_separation(Gc, eff["K_n_Pa_per_m"])

    # Dimensionally-clean physical interface stiffness: through-thickness spring
    # K_n = (lambda+2mu)/eps_phys [Pa/m] (modulus over the dimensional layer
    # thickness). This is regime-independent and yields delta_c ~ nm. It is
    # reported ALONGSIDE the asymptotic-formula value, not as a replacement:
    # the formula K_eff = eta^(alpha-1) K_tilde is only dimensionally Pa/m at
    # alpha=0, so its numeric magnitude (and the resulting delta_c) is convention
    # dependent away from the critical regime. See audit finding KEFF-UNITS.
    K_n_phys = (el["lambda"] + 2.0 * el["mu"]) / eps
    delta_c_phys = P.critical_separation(Gc, K_n_phys)

    return {
        "composition_name": comp_name,
        "eps_m": eps, "L_m": L_MACRO,
        "alpha": round(alpha, 4), "gamma": gamma, "delta": delta,
        "E_sei_Pa": el["E"], "nu_sei": el["nu"],
        "lambda_sei_Pa": el["lambda"], "mu_sei_Pa": el["mu"],
        "D_sei_m2s": D_sei, "beta_sei_m3_per_mol": beta_sei,
        "K_n_eff_Pa_per_m": eff["K_n_Pa_per_m"],
        "K_t_eff_Pa_per_m": eff["K_t_Pa_per_m"],
        "Kn_over_Kt": eff["Kn_over_Kt"],
        "beta_eff_m3_per_mol": eff["beta_eff_m3_per_mol"],
        "R_eff_s_per_m": eff["R_eff_s_per_m"],
        "D_eff_m2s": eff["D_eff_m2s"],
        "g_v_sei_J_per_m3": g_v_sei,
        "Gc_J_per_m2": Gc, "delta_c_m": delta_c,
        "K_n_physical_Pa_per_m": K_n_phys, "delta_c_physical_m": delta_c_phys,
        "mechanical_regime": reg["mechanical_regime"],
        "diffusion_regime": reg["diffusion_regime"],
        "well_posed": reg["well_posed"],
    }


def compute_all_parameters(df: pd.DataFrame) -> pd.DataFrame:
    """Effective parameters for all named compositions, as a DataFrame."""
    return pd.DataFrame([compute_effective_parameters(n, df) for n in SEI_COMPOSITIONS])


def sweep_composition(df: pd.DataFrame, n_points: int = 41) -> pd.DataFrame:
    """
    Sweep the LiF fraction (organic fixed at 0.2, Li2CO3 = balance) and tabulate
    Gc for eps in {10, 30, 50} nm. Drives Figure 4.
    """
    rows = []
    f_org = 0.2
    for f_LiF in np.linspace(0.0, 0.8, n_points):
        f_Li2CO3 = max(0.0, 1.0 - f_LiF - f_org)
        comp = {"LiF": float(f_LiF), "Li2CO3": float(f_Li2CO3), "organic": f_org}
        g_v = homogenize_fracture_energy(df, comp)
        row = {"f_LiF": float(f_LiF), "f_Li2CO3": f_Li2CO3, "f_org": f_org,
               "g_v_J_per_m3": g_v}
        for eps_nm in (10, 30, 50):
            row[f"Gc_{eps_nm}nm_J_per_m2"] = (eps_nm * 1e-9) * g_v
        rows.append(row)
    return pd.DataFrame(rows)
