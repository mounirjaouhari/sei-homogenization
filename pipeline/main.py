#!/usr/bin/env python3
"""
Script pilote principal pour sei_mp_integration.

Usage :
    python main.py                 # Exécute tout le pipeline
    python main.py --offline       # Force le mode offline (sans MP API)
    python main.py --only matlab   # Génère uniquement le fichier MATLAB
    python main.py --only latex    # Génère uniquement les tableaux LaTeX
    python main.py --only figures  # Génère uniquement les figures
    python main.py --sweep         # Balaye les compositions (Figure 4)
"""

import argparse
import sys
from pathlib import Path

# Ajouter le répertoire src au path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from config import (
    SEI_COMPOSITIONS, SEI_THICKNESSES, OUTPUT_DIR,
    OUTPUT_DIR_MATLAB, OUTPUT_DIR_LATEX, OUTPUT_DIR_FIGURES,
    OUTPUT_DIR_DATA, validate_config
)
from fetch_mp_data import fetch_all_constituents
from compute_effective import (
    compute_effective_parameters, compute_all_parameters, sweep_composition
)
from generate_matlab import generate_matlab_config
from generate_latex import (
    generate_table1_constituents, generate_table2_effective,
    generate_section5_tex
)
from plot_composition import (
    plot_gc_vs_composition, plot_phase_diagram,
    plot_effective_params_comparison
)


# =====================================================================
#  Pipeline complet
# =====================================================================

def run_pipeline(use_mp: bool = True,
                 generate_matlab: bool = True,
                 generate_latex: bool = True,
                 generate_figures: bool = True,
                 do_sweep: bool = False):
    """Exécute le pipeline complet de génération."""

    print("=" * 70)
    print("  SEI Materials Project Integration — Pipeline")
    print("=" * 70)

    # 1. Validation de la configuration
    print("\n[1/6] Validation de la configuration...")
    validate_config()
    print("  Configuration valide.")

    # 2. Récupération des données
    print("\n[2/6] Récupération des propriétés des constituants SEI...")
    constituents_df = fetch_all_constituents(use_mp=use_mp)

    # Sauvegarde du DataFrame
    csv_path = OUTPUT_DIR_DATA / "constituents_properties.csv"
    constituents_df.to_csv(csv_path, index=False)
    print(f"  Données sauvegardées : {csv_path}")

    # 3. Calcul des paramètres effectifs
    print("\n[3/6] Calcul des paramètres effectifs (Theorems 3.1 et 4.2)...")
    params_list = []
    for comp_name in SEI_COMPOSITIONS.keys():
        params = compute_effective_parameters(comp_name, constituents_df)
        params_list.append(params)

    # Sauvegarde
    params_df = compute_all_parameters(constituents_df)
    csv_params = OUTPUT_DIR_DATA / "effective_parameters.csv"
    params_df.to_csv(csv_params, index=False)
    print(f"  Paramètres sauvegardés : {csv_params}")

    # 4. Génération du fichier MATLAB
    if generate_matlab:
        print("\n[4/6] Génération du fichier config_parameters_mp.m...")
        matlab_path = generate_matlab_config(params_list)
        print(f"  Fichier MATLAB : {matlab_path}")

    # 5. Génération des tableaux LaTeX
    if generate_latex:
        print("\n[5/6] Génération des tableaux LaTeX...")
        table1 = generate_table1_constituents(constituents_df)
        table2 = generate_table2_effective(params_list)
        section5 = generate_section5_tex(constituents_df, params_list)
        print(f"  Table 1 : {table1}")
        print(f"  Table 2 : {table2}")
        print(f"  Section 5 : {section5}")

    # 6. Génération des figures
    if generate_figures:
        print("\n[6/6] Génération des figures...")
        fig5 = plot_phase_diagram()
        fig6 = plot_effective_params_comparison(params_list)

        if do_sweep:
            print("  Balayage de composition pour Figure 4...")
            sweep_df = sweep_composition(constituents_df)
            sweep_df.to_csv(OUTPUT_DIR_DATA / "composition_sweep.csv", index=False)
            fig4 = plot_gc_vs_composition(sweep_df,
                                          eps_values=[10, 30, 50])
            print(f"  Figure 4 : {fig4}")

    # Résumé
    print("\n" + "=" * 70)
    print("  Pipeline terminé avec succès.")
    print(f"  Répertoires de sortie :")
    print(f"    - MATLAB : {OUTPUT_DIR_MATLAB}")
    print(f"    - LaTeX  : {OUTPUT_DIR_LATEX}")
    print(f"    - Figures: {OUTPUT_DIR_FIGURES}")
    print(f"    - Data   : {OUTPUT_DIR_DATA}")
    print("=" * 70)


# =====================================================================
#  CLI
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="SEI Materials Project Integration"
    )
    parser.add_argument('--offline', action='store_true',
                        help='Force le mode offline (sans MP API)')
    parser.add_argument('--only', choices=['matlab', 'latex', 'figures'],
                        help='Génère uniquement un type de fichier')
    parser.add_argument('--sweep', action='store_true',
                        help='Effectue le balayage de composition (Figure 4)')
    args = parser.parse_args()

    use_mp = not args.offline

    if args.only:
        gen_matlab = (args.only == 'matlab')
        gen_latex = (args.only == 'latex')
        gen_figures = (args.only == 'figures')
    else:
        gen_matlab = gen_latex = gen_figures = True

    run_pipeline(
        use_mp=use_mp,
        generate_matlab=gen_matlab,
        generate_latex=gen_latex,
        generate_figures=gen_figures,
        do_sweep=args.sweep
    )


if __name__ == "__main__":
    main()
