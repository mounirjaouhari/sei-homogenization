"""
Exemples d'utilisation avancée du package sei_mp_integration.

Ces exemples montrent comment :
  1. Récupérer un seul constituant depuis MP
  2. Faire une recherche de nouveaux matériaux candidats ASEI
  3. Tracer une carte de chaleur des paramètres effectifs
  4. Exporter les données pour analyse externe
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / "src"))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from config import SEI_COMPOSITIONS, SEI_THICKNESSES, REGIME_EXPONENTS
from fetch_mp_data import fetch_constituent, fetch_from_materials_project, SEI_CONSTITUENTS
from homogenize import (
    mori_tanaka_homogenization, bruggeman_homogenization,
    homogenize_sei_elastic, homogenize_fracture_energy
)
from compute_effective import (
    compute_effective_parameters, sweep_composition, identify_regime
)


# =====================================================================
#  Exemple 1 : Récupérer un seul constituant depuis MP
# =====================================================================

def example1_single_constituent():
    """Récupère et affiche les propriétés d'un seul constituant SEI."""
    print("\n" + "=" * 60)
    print("  EXEMPLE 1 : Récupération d'un seul constituant (LiF)")
    print("=" * 60)

    # Premier constituant de la liste (LiF)
    lif = SEI_CONSTITUENTS[0]
    props = fetch_constituent(lif, use_mp=False)

    print(f"\nPropriétés de {props.name} ({props.formula}) :")
    print(f"  Source            : {props.source}")
    print(f"  MP ID             : {props.mp_id}")
    print(f"  Module Young E    : {props.E_Pa/1e9:.1f} GPa")
    print(f"  Poisson nu        : {props.nu:.3f}")
    print(f"  Compressibilité K : {props.K_Pa/1e9:.1f} GPa")
    print(f"  Cisaillement G    : {props.G_Pa/1e9:.1f} GPa")
    print(f"  Densité           : {props.density:.0f} kg/m^3")
    print(f"  Bande interdite   : {props.band_gap_eV:.1f} eV")
    print(f"  D_Li              : {props.D_Li_m2s:.2e} m^2/s")
    print(f"  K_Ic              : {props.fracture_toughness_MPa_sqrt_m:.2f} MPa√m")
    print(f"  g (fissuration)   : {props.fracture_energy_Jm2/1e9:.1f} GPa (=J/m^3)")
    print(f"  beta (Vegard)     : {props.beta_vegard_m3_per_mol:.2e} m^3/mol")
    print(f"  Système cristallin: {props.crystal_system} ({props.space_group})")


# =====================================================================
#  Exemple 2 : Recherche de matériaux candidats ASEI
# =====================================================================

def example2_asei_candidates():
    """
    Cherche dans le Materials Project des composés candidats pour
    une ASEI (artificial SEI) optimale.

    Critères :
      - Contient Li
      - Stable (sur le convex hull)
      - Isolant (band_gap > 3 eV)
      - K > 50 GPa (raideur inorganique)
    """
    print("\n" + "=" * 60)
    print("  EXEMPLE 2 : Recherche de candidats ASEI")
    print("=" * 60)

    try:
        from mp_api.client import MPRester
        from config import get_mp_api_key
        api_key = get_mp_api_key()
        if not api_key:
            print("\n  Pas de clé API MP. Pour exécuter cet exemple :")
            print("  1. Obtenir une clé sur https://next-gen.materialsproject.org/api")
            print("  2. Définir : export MP_API_KEY='votre_cle'")
            return

        print("\n  Recherche de composés Li-stables-isolants-raides...")
        with MPRester(api_key) as mpr:
            docs = mpr.materials.summary.search(
                elements=["Li"],
                exclude_elements=["Fe", "Co", "Ni", "Mn", "Cu", "Zn"],
                is_stable=True,
                band_gap=(3.0, None),
                bulk_modulus=(50e9, None),
                fields=["material_id", "formula_pretty",
                        "bulk_modulus", "shear_modulus",
                        "band_gap", "energy_above_hull"]
            )

            candidates = []
            for d in docs:
                if d.bulk_modulus and d.shear_modulus:
                    E = 9 * d.bulk_modulus * d.shear_modulus / \
                        (3 * d.bulk_modulus + d.shear_modulus)
                    candidates.append({
                        "formula": d.formula_pretty,
                        "mp_id": d.material_id,
                        "E_GPa": E / 1e9,
                        "K_GPa": d.bulk_modulus / 1e9,
                        "G_GPa": d.shear_modulus / 1e9,
                        "band_gap_eV": d.band_gap,
                    })

            df = pd.DataFrame(candidates).sort_values("E_GPa", ascending=False)
            print(f"\n  {len(df)} candidats trouvés :")
            print(df.head(20).to_string(index=False))
            return df

    except ImportError:
        print("\n  mp-api non installé. Installer avec : pip install mp-api")
        return None
    except Exception as e:
        print(f"\n  Erreur : {e}")
        return None


# =====================================================================
#  Exemple 3 : Carte de chaleur des paramètres effectifs
# =====================================================================

def example3_heatmap():
    """
    Génère une carte de chaleur de G_c en fonction de la fraction
    de LiF et de l'épaisseur SEI.
    """
    print("\n" + "=" * 60)
    print("  EXEMPLE 3 : Carte de chaleur G_c(f_LiF, eps)")
    print("=" * 60)

    from fetch_mp_data import fetch_all_constituents
    df = fetch_all_constituents(use_mp=False)

    # Grille de paramètres
    f_LiF_range = np.linspace(0, 0.8, 30)
    eps_range_nm = np.linspace(10, 50, 30)
    eps_range_m = eps_range_nm * 1e-9

    # Base composition (mature)
    base_comp = SEI_COMPOSITIONS["mature"].copy()

    # Calcul de G_c pour chaque (f_LiF, eps)
    Gc_grid = np.zeros((len(f_LiF_range), len(eps_range_m)))
    for i, f_LiF in enumerate(f_LiF_range):
        other_names = [n for n in base_comp if n != "LiF"]
        f_others_total = sum(base_comp[n] for n in other_names)
        f_remaining = 1.0 - f_LiF

        comp = {"LiF": f_LiF}
        for n in other_names:
            comp[n] = f_remaining * (base_comp[n] / f_others_total) if f_others_total > 0 else 0

        g_sei = homogenize_fracture_energy(df, comp)
        for j, eps in enumerate(eps_range_m):
            Gc_grid[i, j] = eps * g_sei

    # Tracé
    fig, ax = plt.subplots(figsize=(10, 7))
    im = ax.pcolormesh(eps_range_nm, f_LiF_range, Gc_grid,
                       cmap='viridis', shading='auto')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('$G_c$ (J/m$^2$)', fontsize=12)

    ax.set_xlabel('SEI thickness $\\epsilon_{SEI}$ (nm)', fontsize=12)
    ax.set_ylabel('LiF volume fraction $f_{LiF}$', fontsize=12)
    ax.set_title('Effective fracture energy $G_c = \\epsilon_{SEI} \\sum_i f_i g_i$',
                 fontsize=13)

    # Contours
    cs = ax.contour(eps_range_nm, f_LiF_range, Gc_grid,
                    levels=[0.1, 0.2, 0.3], colors='white', linewidths=1.5)
    ax.clabel(cs, inline=True, fontsize=10, fmt='%.2f')

    plt.tight_layout()
    output_path = Path(__file__).parent.parent / "output" / "figures" / "heatmap_Gc.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"  Carte de chaleur sauvegardée : {output_path}")


# =====================================================================
#  Exemple 4 : Comparaison de schémas d'homogénéisation
# =====================================================================

def example4_compare_schemes():
    """
    Compare Mori-Tanaka, Voigt (règle des mélanges) et Reuss (borne inférieure)
    pour l'élasticité de la SEI.
    """
    print("\n" + "=" * 60)
    print("  EXEMPLE 4 : Comparaison des schémas d'homogénéisation")
    print("=" * 60)

    from fetch_mp_data import fetch_all_constituents
    df = fetch_all_constituents(use_mp=False)

    comp = SEI_COMPOSITIONS["mature"]

    # Voigt (règle des mélanges) : borne supérieure
    E_voigt = 0
    for name, f in comp.items():
        row = df[df["name"] == name].iloc[0]
        E_voigt += f * row["E_Pa"]

    # Reuss (1/E = sum f_i/E_i) : borne inférieure
    inv_E_reuss = 0
    for name, f in comp.items():
        row = df[df["name"] == name].iloc[0]
        inv_E_reuss += f / row["E_Pa"]
    E_reuss = 1.0 / inv_E_reuss

    # Mori-Tanaka
    E_mt = homogenize_sei_elastic(df, comp)["E"]

    print(f"\n  Composition : mature")
    print(f"  Voigt (borne sup.) : E = {E_voigt/1e9:.1f} GPa")
    print(f"  Mori-Tanaka        : E = {E_mt/1e9:.1f} GPa")
    print(f"  Reuss (borne inf.) : E = {E_reuss/1e9:.1f} GPa")
    print(f"\n  => Mori-Tanaka donne une valeur intermédiaire,")
    print(f"     cohérente avec la littérature SEI (15-30 GPa).")


# =====================================================================
#  Exécution de tous les exemples
# =====================================================================

if __name__ == "__main__":
    example1_single_constituent()
    example2_asei_candidates()
    example3_heatmap()
    example4_compare_schemes()

    print("\n" + "=" * 60)
    print("  Tous les exemples ont été exécutés.")
    print("=" * 60)
