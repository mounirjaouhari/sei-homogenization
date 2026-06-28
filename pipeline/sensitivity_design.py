#!/usr/bin/env python3
"""
sensitivity_design.py  --  WP4: robustness (sensitivity) and design (master curve).

(1) SENSITIVITY. The SEI inputs (rupture stress, modulus, diffusivity, fracture
    multiplier) are uncertain by factors of 2-10. We quantify how the headline
    outputs (Gc, delta_c, critical cycling strain, fade@1000) respond, using exact
    log-elasticities  E_X(Y) = d ln Y / d ln X  (computed by perturbation), and a
    tornado over realistic uncertainty ranges. This separates the inputs that
    actually control the predictions from those that do not.

(2) DESIGN MASTER CURVE. The composition-dependent critical cycling strain
        deps_crit(f_org) = sqrt( 2 g_v(f_org) / (lam+2mu)(f_org) )
    is a single dimensionless design curve: a SEI survives cycling iff the cell's
    cyclic strain stays below deps_crit. We sweep the inorganic<->organic balance
    and read off the safe-composition window for graphite / NMC / Si strains.

(3) PREDICTED vs MEASURED. The predicted moduli, Gc and delta_c are compared with
    literature ranges (order-of-magnitude falsifiability).

Run:  python pipeline/sensitivity_design.py
"""

from __future__ import annotations

import sys
from pathlib import Path


def _find_skills():
    for up in Path(__file__).resolve().parents:
        if (up / "Skills" / "sei_physics.py").exists():
            return up / "Skills"
    raise ImportError("Skills/sei_physics.py not found.")


sys.path.insert(0, str(_find_skills()))
import math
import sei_physics as P  # noqa: E402

KAPPA = 15.0
EPS = 30e-9


def props(comp, kappa=KAPPA, sig_scale=1.0, E_scale=1.0):
    """(lam+2mu, g_v) for a composition with optional global scalings of sigma_r, E."""
    fr = P.normalize_fractions(comp)
    # homogenize_elastic uses CONSTITUENTS[k]['E']; E_scale is applied to BOTH the
    # homogenized modulus C1d (below) and to g_v (via E in sigma^2/2E), so the
    # E-elasticity is consistent (-1 for Gc, delta_c and eps_crit).
    el = P.homogenize_elastic(comp)
    C1d = (el["lam"] + 2 * el["mu"]) * E_scale
    g_v = kappa * sum(fr[k] * P.volumetric_fracture_energy(
        P.CONSTITUENTS[k]["sigma_r"] * sig_scale,
        P.CONSTITUENTS[k]["E"] * E_scale) for k in fr)
    return C1d, g_v


def outputs(comp, eps=EPS, kappa=KAPPA, sig_scale=1.0, E_scale=1.0, D_scale=1.0):
    C1d, g_v = props(comp, kappa, sig_scale, E_scale)
    Gc = eps * g_v
    Kn = C1d / eps
    delta_c = math.sqrt(2 * Gc / Kn)
    eps_crit = math.sqrt(2 * g_v / C1d)
    # fade@1000 ~ (eps(N)-eps0); growth k_g ~ D_scale ; eps(N)=sqrt(eps0^2+2 kg tc N)
    eps0, kg_tc0 = 10e-9, ((50e-9) ** 2 - (10e-9) ** 2) / 2000
    epsN = math.sqrt(eps0 ** 2 + 2 * kg_tc0 * D_scale * 1000)
    fade1000 = 1.0 * 5.0e4 * (epsN - eps0) * 96485.0 / (372.0 * 3.6)
    return dict(Gc=Gc, delta_c=delta_c, eps_crit=eps_crit, fade1000=fade1000)


def elasticity(comp, out_key, var, h=1e-3):
    """d ln(output)/d ln(input) by central difference at the nominal point."""
    kw = dict(eps=EPS, kappa=KAPPA, sig_scale=1.0, E_scale=1.0, D_scale=1.0)
    base = {"eps": EPS, "kappa": KAPPA, "sig_scale": 1.0, "E_scale": 1.0, "D_scale": 1.0}
    x = base[var]
    kw[var] = x * (1 + h); yp = outputs(comp, **kw)[out_key]
    kw[var] = x * (1 - h); ym = outputs(comp, **kw)[out_key]
    return (math.log(yp) - math.log(ym)) / (2 * h)


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    comp = P.COMPOSITIONS["mature"]["fractions"]
    print("=" * 72)
    print("  WP4  Robustness (sensitivity) and design (master curve)")
    print("=" * 72)

    # ---- (1) sensitivity: log-elasticities ----
    inputs = [("eps", "SEI thickness"), ("kappa", "fracture multiplier"),
              ("sig_scale", "rupture stress"), ("E_scale", "modulus"), ("D_scale", "growth diffusivity")]
    outs = [("Gc", "Gc"), ("delta_c", "delta_c"), ("eps_crit", "crit. strain"), ("fade1000", "fade@1000")]
    print("\n  Log-elasticities  d ln(Y)/d ln(X)  (mature SEI):")
    print("    {:<20}".format("input \\ output") + "".join("{:>13}".format(o[1]) for o in outs))
    for var, label in inputs:
        row = "    {:<20}".format(label)
        for ok, _ in outs:
            row += "{:>13.2f}".format(elasticity(comp, ok, var))
        print(row)
    print("  => Gc is most sensitive to the rupture stress (elasticity 2, quadratic);")
    print("     the critical strain is independent of thickness; fade scales as sqrt(D).")

    # ---- (2) design master curve: deps_crit vs organic fraction ----
    fracs, dcrit = [], []
    for i in range(1, 20):
        f = i / 20.0
        c = {"Li2CO3": 1 - f, "organic": f}
        _, g_v = props(c); C1d, _ = props(c)
        fracs.append(f); dcrit.append(100 * math.sqrt(2 * g_v / C1d))
    print("\n  Design master curve  deps_crit(f_org)  [%]:")
    for f, d in list(zip(fracs, dcrit))[::4]:
        print(f"    f_org={f:.2f} : crit strain = {d:.2f}%")
    print("  => crack resistance increases with organic fraction (more compliant);")
    print("     trade-off: inorganic passivates (low growth) but cracks earlier ->")
    print("     supports the gradient ASEI (inorganic inner / organic outer).")

    # ---- (3) predicted vs measured ----
    print("\n  Predicted vs measured (order-of-magnitude falsifiability):")
    tab = [
        ("SEI modulus E (GPa)", "2-50 (comp.-dep.)", "0.5-10 (organic, AFM) ... 50-90 (inorganic)"),
        ("Gc (J/m^2)", "0.05-0.32", "0.1-2 (reported SEI fracture energy)"),
        ("delta_c (nm)", "0.4-1.0", "sub-nm to few nm (cohesive zone)"),
        ("crit. strain (%)", "2.0-3.6", "graphite expands ~10% -> cracks (consistent)"),
    ]
    print("    {:<22}{:<14}{}".format("quantity", "predicted", "measured / literature"))
    for q, p, m in tab:
        print("    {:<22}{:<14}{}".format(q, p, m))

    # ---- figures ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        figdir = Path(__file__).resolve().parents[1] / "figures"
        figdir.mkdir(parents=True, exist_ok=True)
        # tornado for fade@1000 and crit strain over realistic uncertainty ranges
        ranges = {"sig_scale": ("rupture stress ×/÷2", 2.0), "E_scale": ("modulus ×/÷2", 2.0),
                  "D_scale": ("growth diff. ×/÷10", 10.0), "kappa": ("fracture mult. ×/÷3", 3.0)}
        fig, axs = plt.subplots(1, 2, figsize=(11, 4))
        for ax, ok, title in [(axs[0], "eps_crit", "critical strain"), (axs[1], "fade1000", "fade @ 1000 cyc")]:
            base = outputs(comp)[ok]
            labels, los, his = [], [], []
            for var, (lab, fac) in ranges.items():
                kw = dict(eps=EPS, kappa=KAPPA, sig_scale=1.0, E_scale=1.0, D_scale=1.0)
                b = {"eps": EPS, "kappa": KAPPA, "sig_scale": 1.0, "E_scale": 1.0, "D_scale": 1.0}
                kw[var] = b[var] * fac; hi = outputs(comp, **kw)[ok]
                kw[var] = b[var] / fac; lo = outputs(comp, **kw)[ok]
                labels.append(lab); los.append(100 * (min(lo, hi) / base - 1)); his.append(100 * (max(lo, hi) / base - 1))
            y = range(len(labels))
            for yi, lo, hi in zip(y, los, his):
                ax.barh(yi, hi - lo, left=lo, color="tab:blue", alpha=0.7)
            ax.set_yticks(list(y)); ax.set_yticklabels(labels, fontsize=8)
            ax.axvline(0, color="k", lw=0.8); ax.set_xlabel("% change in " + title)
            ax.set_title(title)
        fig.tight_layout(); fig.savefig(figdir / "fig_sensitivity.png", dpi=200); plt.close(fig)
        # master curve
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot([100 * f for f in fracs], dcrit, "b-", lw=2.5)
        for strain, name, col in [(2, "NMC ~2%", "green"), (10, "graphite ~10%", "orange")]:
            ax.axhline(strain, ls="--", color=col, lw=1.2, label=f"{name} cycling strain")
        ax.set_xlabel("organic fraction in SEI (%)"); ax.set_ylabel("critical cycling strain Δε_crit (%)")
        ax.fill_between([100 * f for f in fracs], dcrit, 0, color="tab:blue", alpha=0.08)
        ax.set_title("Design master curve: SEI crack resistance vs composition")
        ax.grid(alpha=0.3); ax.legend(loc="center right"); ax.set_ylim(0, 11)
        ax.annotate("graphite always cracks\n(any composition)", (30, 10.2), fontsize=8, color="darkorange")
        ax.annotate("safe for NMC if organic > ~12%", (40, 2.3), fontsize=8, color="darkgreen")
        fig.tight_layout(); fig.savefig(figdir / "fig_master_curve.png", dpi=200); plt.close(fig)
        print(f"\n  Figures written: {figdir/'fig_sensitivity.png'} , {figdir/'fig_master_curve.png'}")
    except Exception as e:  # pragma: no cover
        print(f"  (figures skipped: {e})")
    print("=" * 72)


if __name__ == "__main__":
    main()
