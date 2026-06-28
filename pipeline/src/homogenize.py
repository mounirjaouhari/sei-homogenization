"""
homogenize.py — homogenization of SEI constituents to effective properties.

Delegates the physics to the verified Skills/sei_physics (Mori-Tanaka elastic,
Bruggeman diffusion). The fracture-energy homogenization returns the VOLUMETRIC
density g_v_sei = sum_i f_i g_v^(i) [J/m^3]; multiply by eps_SEI to get Gc [J/m^2].
"""

from __future__ import annotations

from typing import Dict

import pandas as pd

from _engine import P


def _row(df: pd.DataFrame, name: str):
    sub = df[df["name"] == name]
    if sub.empty:
        raise KeyError(f"Constituent '{name}' not found in the constituents table.")
    return sub.iloc[0]


def mori_tanaka_homogenization(E_m: float, nu_m: float, E_i: float, f_i: float) -> float:
    """Mori-Tanaka effective Young's modulus (article Eq. MT)."""
    return P.mori_tanaka_modulus(E_m, nu_m, E_i, f_i)


def bruggeman_homogenization(D_org: float, D_inorg: float, f_inorg: float) -> float:
    """Bruggeman effective diffusivity."""
    return P.bruggeman_diffusivity(D_org, D_inorg, f_inorg)


def homogenize_sei_elastic(df: pd.DataFrame, comp: Dict[str, float]) -> Dict[str, float]:
    """
    Effective (E, nu) of the SEI by Mori-Tanaka (inorganic inclusions in organic
    matrix) + volume-averaged Poisson ratio. comp: {constituent: volume fraction}.
    """
    total = sum(comp.values())
    comp = {k: v / total for k, v in comp.items()}
    organic_keys = [k for k in comp if df.loc[df["name"] == k, "category"].iloc[0] == "organic"]
    f_org = sum(comp[k] for k in organic_keys)
    f_inorg = 1.0 - f_org

    E_m = _row(df, organic_keys[0])["E_Pa"] if organic_keys else _row(df, "organic")["E_Pa"]
    nu_m = _row(df, organic_keys[0])["nu"] if organic_keys else _row(df, "organic")["nu"]

    inorg_keys = [k for k in comp if k not in organic_keys]
    if f_inorg > 0 and inorg_keys:
        E_i = sum(comp[k] * _row(df, k)["E_Pa"] for k in inorg_keys) / f_inorg
        E_eff = mori_tanaka_homogenization(E_m, nu_m, E_i, f_inorg)
    else:
        E_eff = E_m
    nu_eff = sum(comp[k] * _row(df, k)["nu"] for k in comp)
    lam, mu, K = P.lame_from_E_nu(E_eff, nu_eff)
    return {"E": E_eff, "nu": nu_eff, "lambda": lam, "mu": mu, "K": K}


def homogenize_sei_diffusivity(df: pd.DataFrame, comp: Dict[str, float]) -> float:
    """Effective Li diffusivity [m^2/s] by Bruggeman."""
    total = sum(comp.values())
    comp = {k: v / total for k, v in comp.items()}
    organic_keys = [k for k in comp if df.loc[df["name"] == k, "category"].iloc[0] == "organic"]
    f_org = sum(comp[k] for k in organic_keys)
    f_inorg = 1.0 - f_org
    D_org = _row(df, organic_keys[0])["D_Li_m2s"] if organic_keys else _row(df, "organic")["D_Li_m2s"]
    inorg_keys = [k for k in comp if k not in organic_keys]
    if f_inorg > 0 and inorg_keys:
        D_inorg = sum(comp[k] * _row(df, k)["D_Li_m2s"] for k in inorg_keys) / f_inorg
    else:
        D_inorg = D_org
    return bruggeman_homogenization(D_org, D_inorg, f_inorg)


def homogenize_sei_beta(df: pd.DataFrame, comp: Dict[str, float]) -> float:
    """Volume-averaged Vegard coefficient [m^3/mol]."""
    total = sum(comp.values())
    return sum((v / total) * _row(df, k)["beta_vegard_m3_per_mol"] for k, v in comp.items())


def homogenize_fracture_energy(df: pd.DataFrame, comp: Dict[str, float]) -> float:
    """
    Volumetric fracture-energy density g_v_sei = sum_i f_i g_v^(i)  [J/m^3].
    (df['fracture_energy_Jm2'] holds g_v in J/m^3, see fetch_mp_data.)
    Gc = eps_SEI * g_v_sei.
    """
    total = sum(comp.values())
    return sum((v / total) * _row(df, k)["fracture_energy_Jm2"] for k, v in comp.items())
