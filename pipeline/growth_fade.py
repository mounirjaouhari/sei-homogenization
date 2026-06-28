#!/usr/bin/env python3
"""
growth_fade.py  --  WP3: growing SEI -> calibration-free capacity fade.

Two-scale temporal homogenization (slow time tau = eps*t, after Oskay-Fish): the
SEI thickness grows by diffusion-limited (parabolic) kinetics, and the effective
interface parameters evolve with it. The capacity fade is then PREDICTED from
microscopic SEI properties -- it is not fitted to any capacity-vs-cycle curve.

Model
-----
Growth (diffusion-limited through the SEI):   eps d(eps)/dt = k_g
   =>  eps(N) = sqrt(eps0^2 + 2 k_g t_cycle N)        (parabolic, sqrt-N).

Evolving effective parameters (per the validated through-thickness springs):
   K_n(N) = (lam+2mu)/eps(N),   R_eff(N)=eps(N)/D,   Gc(N)=eps(N)*g_v,
   delta_c(N)=sqrt(2 Gc/K_n) = eps(N)*sqrt(2 g_v/(lam+2mu)).

Chemical fade (Li irreversibly built into the growing SEI):
   dQ/Q0 (N) = a_s c_Li,SEI [eps(N)-eps0] F / q_spec        ~ sqrt(N)
   (a_s = electrode specific surface area; c_Li,SEI = Li density in the SEI;
    q_spec = electrode specific capacity).  This is the universal sqrt-t fade.

Mechanical acceleration (composition-dependent threshold):
   the cyclic interfacial energy release rate G_cyc = (1/2)(lam+2mu) deps_v^2 eps
   reaches Gc = g_v eps  iff   deps_v > deps_crit = sqrt(2 g_v/(lam+2mu)),
   INDEPENDENT of eps. Above threshold the SEI cracks and reforms each cycle,
   enhancing the growth (k_g -> k_g (1+chi)); below it, only the chemical sqrt-N
   fade occurs. deps_crit is set by composition through g_v and (lam+2mu).

Run:  python pipeline/growth_fade.py
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

F = 96485.0                      # C/mol

# ---- microscopic, measurable inputs (NOT fitted to capacity data) ----
A_S = 1.0                        # electrode specific surface area [m^2/g] (graphite)
C_LI_SEI = 5.0e4                 # Li density in the SEI solid [mol/m^3]
Q_SPEC = 372.0 * 3.6             # graphite specific capacity 372 mAh/g -> C/g
T_CYCLE = 3600.0                 # s (1 h, ~C/2)
EPS0 = 10e-9                     # initial SEI thickness [m]
# parabolic growth constant calibrated to a MEASURED SEI-thickness trajectory
# (eps: 10 -> 50 nm over 1000 cycles), an independent microscopic observable:
K_G_TCYCLE = ((50e-9) ** 2 - EPS0 ** 2) / (2 * 1000)   # m^2 / cycle


def composition_props(name):
    """(lam+2mu, g_v) for a named SEI composition from the verified physics."""
    comp = P.COMPOSITIONS[name]["fractions"]
    el = P.homogenize_elastic(comp)
    C1d = el["lam"] + 2 * el["mu"]
    # volumetric fracture-energy density (Option B, kappa=15), J/m^3
    kappa = 15.0
    g_v = kappa * sum(P.normalize_fractions(comp)[k] * P.volumetric_fracture_energy(
        P.CONSTITUENTS[k]["sigma_r"], P.CONSTITUENTS[k]["E"]) for k in P.normalize_fractions(comp))
    return C1d, g_v


def eps_of_N(N, kg_tc):
    return math.sqrt(EPS0 ** 2 + 2 * kg_tc * N)


def fade_curve(name, deps_v, N_max=1500):
    """Predicted capacity retention vs cycle for a composition and cycling strain."""
    C1d, g_v = composition_props(name)
    deps_crit = math.sqrt(2 * g_v / C1d)          # composition-dependent critical strain
    cracks = deps_v > deps_crit
    # crack-and-reform growth enhancement, proportional to the strain overshoot
    chi = 2.0 * max(0.0, deps_v / deps_crit - 1.0)
    kg_tc = K_G_TCYCLE * (1 + chi)
    Ns = list(range(0, N_max + 1, 10))
    retention = []
    for N in Ns:
        eps = eps_of_N(N, kg_tc)
        dQ = A_S * C_LI_SEI * (eps - EPS0) * F / Q_SPEC      # cumulative fade fraction
        retention.append(100.0 * (1.0 - dQ))
    return Ns, retention, deps_crit, cracks


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    print("=" * 70)
    print("  Calibration-free capacity fade from a growing SEI (WP3)")
    print("=" * 70)
    print(f"  Inputs (microscopic, not fitted): a_s={A_S} m2/g, c_Li,SEI={C_LI_SEI:.0e} mol/m3,")
    print(f"  q_spec={Q_SPEC:.0f} C/g, eps:10->50 nm/1000cyc (measured growth).")
    print()
    scenarios = [("aged", 0.03), ("mature", 0.03), ("young", 0.03)]   # 3% cycling strain
    curves = {}
    for name, dev in scenarios:
        Ns, ret, dcrit, cr = fade_curve(name, dev)
        curves[(name, dev)] = (Ns, ret)
        print(f"  {name:7} (cycling strain {100*dev:.0f}%): critical strain "
              f"{100*dcrit:.2f}%  -> {'CRACKS (accelerated)' if cr else 'no cracking (sqrt-N only)'}")
        print(f"            retention @ 500 cyc = {ret[Ns.index(500)]:.1f}% , "
              f"@1000 = {ret[Ns.index(1000)]:.1f}%")
    print()
    print("  => fade ~ sqrt(N) (universal SEI-growth signature), with a composition-")
    print("     dependent cracking threshold that accelerates it. No capacity-curve fit.")

    # figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 5))
        styles = {"aged": ("tab:red", "-"), "mature": ("tab:blue", "-"), "young": ("tab:green", "--")}
        for (name, dev), (Ns, ret) in curves.items():
            c, ls = styles[name]
            ax.plot(Ns, ret, ls, color=c, lw=2, label=f"{name} (Δε_v={100*dev:.0f}%)")
        ax.set_xlabel("cycle number N"); ax.set_ylabel("capacity retention (%)")
        ax.set_title("Predicted capacity fade from growing SEI (calibration-free)")
        ax.grid(alpha=0.3); ax.legend(); ax.set_ylim(70, 100)
        out = Path(__file__).resolve().parents[1] / "figures" / "fig_capacity_fade.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)
        print(f"\n  Figure written: {out}")
    except Exception as e:  # pragma: no cover
        print(f"  (figure skipped: {e})")
    print("=" * 70)


if __name__ == "__main__":
    main()
