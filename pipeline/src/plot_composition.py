"""
plot_composition.py — figures (headless Agg backend).

  - plot_gc_vs_composition : Gc vs LiF fraction for three thicknesses (Figure 4)
  - plot_phase_diagram     : (alpha, gamma) phase diagram with SEI markers
  - plot_effective_params_comparison : bar chart of K_n / R_eff / Gc
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Dict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from config import OUTPUT_DIR_FIGURES, REGIME_EXPONENTS


def plot_gc_vs_composition(sweep_df, eps_values=(10, 30, 50),
                           out_dir: Path = OUTPUT_DIR_FIGURES) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "fig4_Gc_vs_composition.png"
    fig, ax = plt.subplots(figsize=(8, 5.5))
    colors = {10: "tab:orange", 30: "tab:blue", 50: "tab:red"}
    for eps_nm in eps_values:
        col = f"Gc_{eps_nm}nm_J_per_m2"
        ax.plot(sweep_df["f_LiF"], sweep_df[col], lw=2.2,
                color=colors.get(eps_nm), label=fr"$\epsilon_{{SEI}}={eps_nm}$ nm")
    ax.axhspan(0, 0.1, color="grey", alpha=0.12)
    ax.text(0.02, 0.05, r"$G_c<0.1$ J/m$^2$ (brittle)", fontsize=9, color="0.3")
    ax.set_xlabel(r"LiF volume fraction $f_{\mathrm{LiF}}$", fontsize=12)
    ax.set_ylabel(r"$G_c$ (J/m$^2$)", fontsize=12)
    ax.set_title(r"Predicted SEI fracture energy $G_c=\epsilon_{SEI}\sum_i f_i g_{v,i}$")
    ax.set_xlim(0, 0.8); ax.set_ylim(0, 0.6)
    ax.grid(alpha=0.25); ax.legend()
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)
    return path


def plot_phase_diagram(params_list: List[Dict] = None,
                       out_dir: Path = OUTPUT_DIR_FIGURES) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "fig_phase_diagram.png"
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for x in (-1, 0, 1):
        ax.axvline(x, ls="--", color="tab:red", alpha=0.5)
    for y in (1, 2):
        ax.axhline(y, ls="--", color="tab:blue", alpha=0.5)
    labels = {"soft": 1.5, "intermediate": 0.5, "critical": -0.5, "stiff": -1.5}
    for name, xc in labels.items():
        ax.text(xc, 2.6, name, ha="center", color="tab:red", fontsize=9)
    markers = params_list or [
        {"composition_name": n, "alpha": REGIME_EXPONENTS[n]["alpha"],
         "gamma": REGIME_EXPONENTS[n]["gamma"]} for n in REGIME_EXPONENTS]
    for m in markers:
        ax.plot(m["alpha"], m["gamma"], "*", ms=15, color="orange",
                markeredgecolor="k")
        ax.annotate(m["composition_name"], (m["alpha"], m["gamma"]),
                    textcoords="offset points", xytext=(6, 4), fontsize=9)
    ax.set_xlabel(r"Elastic exponent $\alpha$", fontsize=12)
    ax.set_ylabel(r"Diffusion exponent $\gamma$", fontsize=12)
    ax.set_title("Phase diagram of SEI asymptotic regimes")
    ax.set_xlim(-2, 2); ax.set_ylim(-0.5, 3); ax.grid(alpha=0.2)
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)
    return path


def plot_effective_params_comparison(params_list: List[Dict],
                                     out_dir: Path = OUTPUT_DIR_FIGURES) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "fig_effective_params.png"
    names = [p["composition_name"] for p in params_list]
    Kn = [p["K_n_eff_Pa_per_m"] / 1e9 * 1e-6 for p in params_list]  # GPa/um
    Gc = [p["Gc_J_per_m2"] for p in params_list]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    ax1.bar(names, Kn, color="tab:blue", alpha=0.8)
    ax1.set_ylabel(r"$K_n$ (GPa/$\mu$m)"); ax1.set_title("Normal interface stiffness")
    ax1.grid(alpha=0.25, axis="y")
    ax2.bar(names, Gc, color="tab:green", alpha=0.8)
    ax2.set_ylabel(r"$G_c$ (J/m$^2$)"); ax2.set_title("Effective fracture energy")
    ax2.grid(alpha=0.25, axis="y")
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)
    return path
