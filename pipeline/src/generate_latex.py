"""
generate_latex.py — emit LaTeX snippets (Table 1, Table 3, Section-5 numbers)
consistent with the CORRECTED article (volumetric g_v column, Gc 0.05-0.32).
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Dict

import pandas as pd

from config import OUTPUT_DIR_LATEX


def generate_table1_constituents(df: pd.DataFrame,
                                 out_dir: Path = OUTPUT_DIR_LATEX) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "table1_constituents.tex"
    L = []
    a = L.append
    a("% Auto-generated Table 1 (constituents). Volumetric g_v in MJ/m^3.")
    a(r"\begin{tabular}{@{}l c c c c c@{}}")
    a(r"\toprule")
    a(r"Constituent & $E$ (GPa) & $\nu$ & $\sigma_r$ (MPa) & $g_v$ (MJ/m$^3$) & $D_{\mathrm{Li}^+}$ (m$^2$/s) \\")
    a(r"\midrule")
    for _, r in df.iterrows():
        a(f"{_tex_name(r['name'])} & {r['E_Pa']/1e9:.3g} & {r['nu']:.2f} & "
          f"{r['sigma_r_Pa']/1e6:.0f} & {r['fracture_energy_Jm2']/1e6:.1f} & "
          f"${_sci(r['D_Li_m2s'])}$ \\\\")
    a(r"\bottomrule")
    a(r"\end{tabular}")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")
    return path


def generate_table2_effective(params_list: List[Dict],
                              out_dir: Path = OUTPUT_DIR_LATEX) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "table2_effective.tex"
    L = []
    a = L.append
    a("% Auto-generated Table 3 (predicted effective parameters).")
    a(r"\begin{tabular}{@{}l c c c c c c@{}}")
    a(r"\toprule")
    a(r"SEI & $\epsS$ (nm) & Regime & $K_n$ (GPa/$\mu$m) & $R_{\mathrm{eff}}$ (s/m) & "
      r"$\Gc$ (J/m$^2$) & $\delta_c$ (nm) \\")
    a(r"\midrule")
    for p in params_list:
        Kn_GPa_per_um = p["K_n_eff_Pa_per_m"] / 1e9 * 1e-6   # Pa/m -> GPa/um
        a(f"{p['composition_name'].capitalize()} & {p['eps_m']*1e9:.0f} & "
          f"{p['mechanical_regime']} & {Kn_GPa_per_um:.3g} & "
          f"${_sci(p['R_eff_s_per_m'])}$ & {p['Gc_J_per_m2']:.3f} & "
          f"{p['delta_c_m']*1e9:.2f} \\\\")
    a(r"\bottomrule")
    a(r"\end{tabular}")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")
    return path


def generate_section5_tex(df: pd.DataFrame, params_list: List[Dict],
                          out_dir: Path = OUTPUT_DIR_LATEX) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "section5_numbers.tex"
    L = []
    a = L.append
    a("% Auto-generated Section-5 numerical values (corrected conventions).")
    for p in params_list:
        a(f"% {p['composition_name']}: alpha={p['alpha']:.3f} "
          f"({p['mechanical_regime']}), E_SEI={p['E_sei_Pa']/1e9:.1f} GPa, "
          f"Gc={p['Gc_J_per_m2']:.3f} J/m^2, delta_c={p['delta_c_m']*1e9:.2f} nm")
    a("")
    a(r"\newcommand{\GcYoung}{" + f"{params_list[0]['Gc_J_per_m2']:.2f}" + r"}")
    a(r"\newcommand{\GcMature}{" + f"{params_list[1]['Gc_J_per_m2']:.2f}" + r"}")
    a(r"\newcommand{\GcAged}{" + f"{params_list[2]['Gc_J_per_m2']:.2f}" + r"}")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")
    return path


# --------------------------------------------------------------------------
def _tex_name(name: str) -> str:
    return {"Li2CO3": r"Li$_2$CO$_3$", "Li2O": r"Li$_2$O",
            "organic": "Alkylcarb.", "polyolefins": "Polyolefins"}.get(name, name)


def _sci(x: float) -> str:
    if x == 0:
        return "0"
    import math
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    return f"{m:.1f}\\times10^{{{e}}}"
