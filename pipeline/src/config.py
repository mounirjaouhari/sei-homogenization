"""
config.py — configuration of the SEI Materials-Project integration pipeline.

All numerical conventions match the CORRECTED article (main.tex):
  - power laws use the dimensionless base eta = eps/L;
  - the fracture-energy density is g_v = KAPPA * sigma_r^2/(2E) [J/m^3]
    (Option B, process-zone multiplier KAPPA), giving Gc ~ 0.05-0.32 J/m^2.
"""

from __future__ import annotations

import os
from pathlib import Path

# --------------------------------------------------------------------------
#  Convention constants (kept consistent with main.tex / Skills/sei_physics)
# --------------------------------------------------------------------------
L_MACRO = 10e-6          # macroscopic length scale [m] (eta = eps/L)
E_ELECTRODE = 30e9       # reference electrode modulus Etilde [Pa]
KAPPA = 15.0             # process-zone multiplier: g_v = KAPPA * sigma_r^2/(2E)

# --------------------------------------------------------------------------
#  Representative SEI compositions (volume fractions, sum to 1) and thickness
#  Mirrors Skills/sei_physics.COMPOSITIONS so every layer agrees.
# --------------------------------------------------------------------------
SEI_COMPOSITIONS = {
    "young":  {"LiF": 0.20, "Li2CO3": 0.30, "Li2O": 0.05, "organic": 0.45},
    "mature": {"LiF": 0.30, "Li2CO3": 0.40, "Li2O": 0.10, "organic": 0.20},
    "aged":   {"LiF": 0.45, "Li2CO3": 0.35, "Li2O": 0.12, "organic": 0.08},
}

SEI_THICKNESSES = {"young": 10e-9, "mature": 30e-9, "aged": 50e-9}

# Corrected scaling exponents (alpha identified with the dimensionless base eta;
# gamma, delta as in the article). alpha here is indicative; compute_effective
# re-identifies it from the homogenized modulus via identify_regime().
REGIME_EXPONENTS = {
    "young":  {"alpha": +0.26, "gamma": 1.5, "delta": 0.0},
    "mature": {"alpha": +0.07, "gamma": 1.2, "delta": 0.0},
    "aged":   {"alpha": -0.03, "gamma": 0.8, "delta": 0.0},
}

# --------------------------------------------------------------------------
#  Output directories (under programmes/output, non-destructive)
# --------------------------------------------------------------------------
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "output"
OUTPUT_DIR_MATLAB = OUTPUT_DIR / "matlab"
OUTPUT_DIR_LATEX = OUTPUT_DIR / "latex"
OUTPUT_DIR_FIGURES = OUTPUT_DIR / "figures"
OUTPUT_DIR_DATA = OUTPUT_DIR / "data"


def get_mp_api_key():
    """Return the Materials Project API key from the environment, or None."""
    return os.environ.get("MP_API_KEY") or None


def validate_config():
    """Sanity-check the configuration; create output directories."""
    for name, comp in SEI_COMPOSITIONS.items():
        s = sum(comp.values())
        if abs(s - 1.0) > 1e-9:
            raise ValueError(f"Composition '{name}' does not sum to 1 (got {s}).")
        if name not in SEI_THICKNESSES:
            raise ValueError(f"Missing thickness for composition '{name}'.")
        if name not in REGIME_EXPONENTS:
            raise ValueError(f"Missing exponents for composition '{name}'.")
    for d in (OUTPUT_DIR, OUTPUT_DIR_MATLAB, OUTPUT_DIR_LATEX,
              OUTPUT_DIR_FIGURES, OUTPUT_DIR_DATA):
        d.mkdir(parents=True, exist_ok=True)
    return True
