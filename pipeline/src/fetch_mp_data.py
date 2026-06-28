"""
fetch_mp_data.py — constituent properties (offline by default).

Offline mode returns curated literature values (identical to the verified
Skills/sei_physics constituent table), augmented with Materials-Project-style
metadata. If `use_mp=True` and `mp-api` + an `MP_API_KEY` are available, the
elastic moduli / band gap are refreshed from the Materials Project; otherwise the
literature values are used and a note is printed.

The `fracture_energy_Jm2` field stores the VOLUMETRIC fracture-energy density
g_v = KAPPA * sigma_r^2/(2E) in J/m^3 (Option B), consistent with the corrected
article and with Gc = eps_SEI * sum_i f_i g_v^(i).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import List, Optional

import pandas as pd

from _engine import P
from config import KAPPA, get_mp_api_key


@dataclass
class ConstituentProperties:
    name: str
    formula: str
    mp_id: str
    category: str
    source: str
    E_Pa: float
    nu: float
    K_Pa: float
    G_Pa: float
    density: float
    band_gap_eV: float
    D_Li_m2s: float
    sigma_r_Pa: float
    fracture_toughness_MPa_sqrt_m: float
    fracture_energy_Jm2: float      # volumetric g_v [J/m^3] (see module docstring)
    beta_vegard_m3_per_mol: float
    crystal_system: str
    space_group: str


# Curated descriptors (literature / Materials Project). sigma_r [Pa], E/nu/D/beta
# come from the verified Skills/sei_physics.CONSTITUENTS table.
_META = {
    "LiF":     dict(formula="LiF",      mp_id="mp-1137",   category="inorganic",
                    K_Pa=35.6e9, G_Pa=26.2e9, density=2640.0, band_gap_eV=13.7,
                    crystal_system="cubic", space_group="Fm-3m"),
    "Li2CO3":  dict(formula="Li2CO3",   mp_id="mp-755798", category="inorganic",
                    K_Pa=39.4e9, G_Pa=23.6e9, density=2110.0, band_gap_eV=5.4,
                    crystal_system="monoclinic", space_group="C2/m"),
    "Li2O":    dict(formula="Li2O",     mp_id="mp-1960",   category="inorganic",
                    K_Pa=44.4e9, G_Pa=33.3e9, density=2010.0, band_gap_eV=7.7,
                    crystal_system="cubic", space_group="Fm-3m"),
    "organic": dict(formula="Li2C2O4",  mp_id="mp-23953",  category="organic",
                    K_Pa=2.2e9, G_Pa=0.74e9, density=1980.0, band_gap_eV=3.5,
                    crystal_system="monoclinic", space_group="P21/c"),
    "polyolefins": dict(formula="(C2H4)n", mp_id="",        category="organic",
                    K_Pa=2.0e9, G_Pa=0.6e9, density=950.0, band_gap_eV=6.0,
                    crystal_system="amorphous", space_group="-"),
}

# Ordered list of constituents to expose (LiF first).
SEI_CONSTITUENTS: List[str] = ["LiF", "Li2CO3", "Li2O", "organic", "polyolefins"]


def _volumetric_g(sigma_r: float, E: float) -> float:
    """g_v = KAPPA * sigma_r^2/(2E)  [J/m^3]  (Option B)."""
    return KAPPA * P.volumetric_fracture_energy(sigma_r, E)


def fetch_constituent(name: str, use_mp: bool = False) -> ConstituentProperties:
    """Return the properties of one SEI constituent (offline literature by default)."""
    if name not in _META:
        raise KeyError(f"Unknown constituent '{name}'. Known: {SEI_CONSTITUENTS}")
    meta = _META[name]
    # Base physics from the verified single source of truth.
    base = P.CONSTITUENTS.get(name) or P.CONSTITUENTS["organic"]
    E, nu = base["E"], base["nu"]
    sigma_r = base["sigma_r"]
    K_Ic = base["K_Ic"]
    D = base["D"]
    beta = base["beta"]
    source = "literature"

    if use_mp:
        mp = _try_materials_project(meta.get("mp_id", ""))
        if mp is not None:
            E = mp.get("E_Pa", E)
            K = mp.get("K_Pa", meta["K_Pa"])
            G = mp.get("G_Pa", meta["G_Pa"])
            meta = {**meta, "K_Pa": K, "G_Pa": G,
                    "band_gap_eV": mp.get("band_gap_eV", meta["band_gap_eV"])}
            source = "materials_project"

    return ConstituentProperties(
        name=name, formula=meta["formula"], mp_id=meta["mp_id"],
        category=meta["category"], source=source,
        E_Pa=E, nu=nu, K_Pa=meta["K_Pa"], G_Pa=meta["G_Pa"],
        density=meta["density"], band_gap_eV=meta["band_gap_eV"], D_Li_m2s=D,
        sigma_r_Pa=sigma_r, fracture_toughness_MPa_sqrt_m=K_Ic / 1e6,
        fracture_energy_Jm2=_volumetric_g(sigma_r, E),
        beta_vegard_m3_per_mol=beta,
        crystal_system=meta["crystal_system"], space_group=meta["space_group"],
    )


def fetch_from_materials_project(mp_id: str) -> Optional[dict]:
    """Public wrapper around the (optional) Materials Project lookup."""
    return _try_materials_project(mp_id)


def _try_materials_project(mp_id: str) -> Optional[dict]:
    """Try to fetch elastic data from MP; return None if unavailable."""
    if not mp_id:
        return None
    api_key = get_mp_api_key()
    if not api_key:
        return None
    try:  # pragma: no cover - network/optional
        from mp_api.client import MPRester
        with MPRester(api_key) as mpr:
            doc = mpr.materials.summary.search(
                material_ids=[mp_id],
                fields=["bulk_modulus", "shear_modulus", "band_gap"])[0]
            K, G = doc.bulk_modulus, doc.shear_modulus
            if K and G:
                E = 9 * K * G / (3 * K + G)
                return dict(E_Pa=E, K_Pa=K, G_Pa=G, band_gap_eV=doc.band_gap)
    except Exception as e:  # pragma: no cover
        print(f"  [MP] lookup failed for {mp_id}: {e}")
    return None


def fetch_all_constituents(use_mp: bool = False) -> pd.DataFrame:
    """Return a DataFrame of all SEI constituents."""
    if use_mp and not get_mp_api_key():
        print("  [fetch] No MP_API_KEY -> using offline literature values.")
        use_mp = False
    rows = [asdict(fetch_constituent(n, use_mp=use_mp)) for n in SEI_CONSTITUENTS]
    return pd.DataFrame(rows)
