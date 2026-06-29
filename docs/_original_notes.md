# SEI Asymptotic Validation — MATLAB Program

Validation numérique des conditions de saut chémo-mécaniques dérivées par
expansion asymptotique pour la SEI (Solid Electrolyte Interphase) dans les
batteries lithium-ion, tel que décrit dans l'article :

> *Asymptotic Derivation of Chemo-Mechanical Jump Conditions for the Solid
> Electrolyte Interphase in Lithium-Ion Batteries* — Section 6.

## Objectif

Le programme compare deux modèles pour une configuration planaire
électrode/SEI/électrolyte soumise à un couplage chémo-mécanique :

1. **Modèle complet (FULL)** — la SEI est maillée explicitement avec
   suffisamment d'éléments dans l'épaisseur (≥ 10 éléments quadratiques).
2. **Modèle effectif (EFF)** — la SEI est remplacée par une interface
   d'épaisseur nulle portant les conditions de saut dérivées
   (Theorem 3.1 de l'article).

L'objectif est de valider la convergence théorique en O(ε) de l'erreur
relative entre les deux modèles pour les quatre régimes asymptotiques
(souple, intermédiaire, critique, raide).

## Architecture

```
sei_validation_matlab/
├── README.md                        Cette documentation
├── main_validation.m                Script pilote principal
│
├── config/
│   └── config_parameters.m          Paramètres matériaux et régimes
│
├── mesh/
│   ├── build_trilayer_mesh.m        Maillage 1D du modèle FULL
│   └── build_eff_mesh.m             Maillage 1D du modèle EFF
│
├── assembly/
│   ├── elastic_element_stiff.m      Matrice de raideur élémentaire (élastique)
│   ├── diffusion_element_stiff.m    Matrice de raideur élémentaire (diffusion)
│   ├── assemble_full_model.m        Assemblage FULL (élastique + diffusion)
│   ├── assemble_eff_model.m         Assemblage EFF avec conditions de saut
│   └── interface_element.m          Élément interface chémo-mécanique
│
├── solve/
│   ├── apply_dirichlet_bc.m         Application des conditions de Dirichlet
│   ├── solve_steady_state.m         Solveur stationnaire
│   └── solve_transient_cyclic.m     Solveur transitoire (Euler implicite)
│
├── cohesive/
│   ├── chemo_cohesive_law.m         Loi cohesive chémo-mécanique (régime critique)
│   └── compute_damage_field.m       Calcul du champ d'endommagement
│
├── post/
│   ├── compute_l2_error.m           Erreur L2 relative
│   ├── plot_convergence_study.m     Courbes de convergence (Figure 5 de l'article)
│   ├── plot_field_comparison.m      Comparaison de champs (Figure 1 de l'article)
│   └── plot_cohesive_curves.m       Lois cohesives (Figure 3 de l'article)
│
├── tests/
│   ├── test_case1_steady.m          Cas 1 : trilayer plan, régime stationnaire
│   ├── test_case2_transient.m       Cas 2 : chargement cyclique C/2
│   ├── test_case3_fracture.m        Cas 3 : endommagement et propagation de fissure
│   └── test_case4_convergence.m     Cas 4 : étude de convergence tous régimes
│
└── utils/
    ├── effective_params.m           Calcul de K_eff, beta_eff, R_eff, G_c
    ├── acoustic_tensor.m            Tenseur acoustique pour matériau isotrope
    └── safe_spmatrix.m              Construction robuste de matrices creuses
```

## Modélisation

Le programme utilise un modèle 1D à travers l'épaisseur (direction z),
suffisant pour capturer la physique essentielle de la SEI :

- Élasticité linéaire 1D avec couplage chémo-mécanique (eigenstrain β·c)
- Diffusion transitoire de Li avec coefficient effectif D
- Élément fini quadratique (3 nœuds par élément)
- Conditions de saut dérivées : K_eff·[[u]] = σ·n - β_eff·<c>,  [[c]] = R_eff·J

## Paramètres

Les paramètres matériaux proviennent de la Table 1 de l'article
(constituants SEI : LiF, Li₂CO₃, Li₂O, alkylcarbonates, polyoléfines).

## Utilisation

```matlab
>> cd sei_validation_matlab
>> main_validation         % Lance tous les cas de test
```

Pour un test individuel :

```matlab
>> addpath(genpath(pwd))
>> test_case4_convergence  % Étude de convergence uniquement
```

## Sorties

- `results/case1_*.mat`        — Données du cas 1 (steady-state)
- `results/case2_*.mat`        — Données du cas 2 (transitoire)
- `results/case3_*.mat`        — Données du cas 3 (fracture)
- `results/case4_convergence.mat` — Étude de convergence
- `figures/fig1_geometry.png`      — Schéma géométrique
- `figures/fig2_phase_diagram.png` — Diagramme de phases
- `figures/fig3_cohesive.png`      — Lois cohesives
- `figures/fig5_convergence.png`   — Convergence O(ε)
- `figures/fig6_fracture.png`      — Comparaison fissure

## Auteur

J.Mounir, École Normale Supérieur CASABLANCA.
