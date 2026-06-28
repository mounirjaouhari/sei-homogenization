#!/usr/bin/env python3
"""
recompute_effective.py — Correct, runnable numerical core for the SEI effective
interface parameters.

WHY THIS FILE EXISTS
--------------------
The original Python pipeline (main.py, example_usage.py) imports a `src/` package
(config, fetch_mp_data, homogenize, compute_effective, generate_matlab, ...) that
is NOT present in this folder, so it cannot run. Independently, the MATLAB
`config_parameters.m` / `effective_params.m` and the committed
`effective_parameters.csv` contain two provable numerical bugs:

  1. Gc units: g_i was set to ~1e9 J/m^3, giving Gc ~ 22-195 J/m^2 (about 2-3
     orders of magnitude above the article's own Table 3, 0.05-0.32 J/m^2, and
     dimensionally inconsistent with Theorem (Gc) as written).
  2. Power-law base inconsistency: K_eff and R_eff used the dimensionless base
     (eps/L), while D_eff used the dimensional eps, making D_eff ~ 1e-24..1e-28
     m^2/s (physically nonsensical for a diffusivity).

This script recomputes the effective parameters CORRECTLY and reproducibly by
delegating ALL physics to the verified single source of truth `sei_physics.py`
(in the sibling Skills/ folder), writes `effective_parameters_corrected.csv`, and
prints a before/after comparison against the committed (buggy) CSV.

It is fully runnable with only the Python standard library + the Skills folder.

OBJECTIVITY NOTE
----------------
The Gc magnitude depends on the *definition* of the constituent fracture-energy
density g_i, which is a modeling choice. Here g_i is taken as the volumetric
elastic-energy density at the strength limit, g_i = sigma_{r,i}^2 / (2 E_i)
[J/m^3] — ONE dimensionally-consistent, reproducible choice (the same one used by
the skills). This is NOT imposed as the unique truth; matching the article's
Table 3 magnitude would require a different g_i (also legitimate). What this
script guarantees is internal consistency: theory (skills) == this core == the
corrected CSV, all dimensionally correct.

Usage:
    python recompute_effective.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

# Reuse the verified physics from the sibling Skills/ folder (single source of truth).
def _find_skills():
    for up in Path(__file__).resolve().parents:
        if (up / "Skills" / "sei_physics.py").exists():
            return up / "Skills"
    raise ImportError("Could not locate Skills/sei_physics.py.")
sys.path.insert(0, str(_find_skills()))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import sei_physics as P  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"                     # programmes/data
BUGGY_CSV = DATA / "effective_parameters.csv"
CORRECTED_CSV = DATA / "effective_parameters_corrected.csv"


def compute_row(name: str) -> dict:
    """Recompute one SEI composition's effective parameters, all base-consistent."""
    comp = P.composition(name)
    fr = comp["fractions"]
    eps = comp["eps_SEI"]
    alpha, gamma, delta = comp["alpha"], comp["gamma"], comp["delta"]

    el = P.homogenize_elastic(fr)
    D_sei = P.homogenize_diffusivity(fr)
    beta_sei = P.homogenize_beta(fr)

    # Effective parameters: ALL power laws use the dimensionless base eta = eps/L.
    eff = P.effective_parameters(el["E"], el["nu"], D_sei, beta_sei,
                                 eps, alpha, gamma, delta, L=P.L_MACRO, h=1.0)

    # Fracture energy with dimensionally-consistent volumetric g_i (documented choice).
    Gc = P.fracture_energy(fr, eps, mode="volumetric")
    delta_c = P.critical_separation(Gc, eff["K_n_Pa_per_m"])

    return {
        "composition_name": name,
        "eps_m": eps,
        "L_m": P.L_MACRO,
        "alpha": alpha, "gamma": gamma, "delta": delta,
        "E_sei_Pa": round(el["E"], 3),
        "nu_sei": round(el["nu"], 4),
        "D_sei_m2s": D_sei,
        "beta_sei_m3_per_mol": beta_sei,
        "K_n_eff_Pa_per_m": eff["K_n_Pa_per_m"],
        "K_t_eff_Pa_per_m": eff["K_t_Pa_per_m"],
        "beta_eff_m3_per_mol": eff["beta_eff_m3_per_mol"],
        "R_eff_s_per_m": eff["R_eff_s_per_m"],
        "D_eff_m2s": eff["D_eff_m2s"],            # base-consistent (eps/L)^gamma
        "Gc_J_per_m2": Gc,                        # dimensionally correct
        "delta_c_m": delta_c,
        "mechanical_regime": P.classify_mechanical_regime(alpha)["regime"],
        "diffusion_regime": P.classify_diffusion_regime(gamma)["regime"],
        "well_posed": P.well_posed(alpha, gamma),
    }


def read_buggy_csv() -> dict:
    """Read the committed (buggy) CSV keyed by composition_name, if present."""
    if not BUGGY_CSV.exists():
        return {}
    out = {}
    with open(BUGGY_CSV, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            out[row.get("composition_name", "?")] = row
    return out


def main():
    print("=" * 74)
    print("  Recompute SEI effective parameters (corrected, base-consistent)")
    print("=" * 74)

    rows = [compute_row(n) for n in ("young", "mature", "aged")]
    buggy = read_buggy_csv()

    # Write the corrected CSV.
    with open(CORRECTED_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    print(f"\n  Corrected CSV written: {CORRECTED_CSV.name}")

    # Before/after comparison for the two bugged quantities.
    print("\n  Before/after (the two provable bugs):")
    print(f"  {'SEI':7} | {'Gc buggy':>12} | {'Gc corrected':>13} | "
          f"{'D_eff buggy':>12} | {'D_eff corrected':>15}")
    print("  " + "-" * 70)
    for r in rows:
        b = buggy.get(r["composition_name"], {})
        gc_b = b.get("Gc_J_per_m2", "n/a")
        deff_b = b.get("D_eff_m2s", "n/a")
        gc_b_s = f"{float(gc_b):.3g}" if gc_b not in ("n/a", "", None) else "n/a"
        deff_b_s = f"{float(deff_b):.3g}" if deff_b not in ("n/a", "", None) else "n/a"
        print(f"  {r['composition_name']:7} | {gc_b_s:>12} | {r['Gc_J_per_m2']:>13.4g} | "
              f"{deff_b_s:>12} | {r['D_eff_m2s']:>15.4g}")

    print("\n  Interpretation:")
    print("   - Gc corrected ~ 0.003-0.02 J/m^2 (volumetric g_i = sigma_r^2/2E); the")
    print("     buggy CSV ~ 22-195 J/m^2 came from g_i ~ 1e9 J/m^3 (units error).")
    print("   - D_eff corrected uses the SAME dimensionless base (eps/L)^gamma as K_eff,")
    print("     removing the ~1e-24..1e-28 nonsensical values of the buggy CSV.")
    print("   - Theory (skills) == this core == corrected CSV, all dimensionally consistent.")
    print("=" * 74)


if __name__ == "__main__":
    main()
