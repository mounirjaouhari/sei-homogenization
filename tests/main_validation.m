%========================================================================
%  MAIN_VALIDATION  Programme principal de validation numérique
%
%  Valide les conditions de saut chémo-mécaniques dérivées par expansion
%  asymptotique pour la SEI dans les batteries lithium-ion.
%
%  Référence : Section 6 de l'article
%    "Asymptotic Derivation of Chemo-Mechanical Jump Conditions for the
%     Solid Electrolyte Interphase in Lithium-Ion Batteries"
%
%  Auteurs : [Author], J.-J. Marigo
%  Date    : 2026
%
%  Architecture :
%    - Cas 1 : trilayer stationnaire (validation de base, 1 régime)
%    - Cas 2 : chargement cyclique transitoire (validation temporelle)
%    - Cas 3 : endommagement et fissuration (régime critique)
%    - Cas 4 : étude de convergence sur les 4 régimes (Figure 5)
%
%  Utilisation :
%    >> cd sei_validation_matlab
%    >> main_validation
%
%  Pour un cas individuel :
%    >> addpath(genpath(pwd))
%    >> P = config_parameters('intermediate');
%    >> test_case1_steady(P);
%========================================================================

function main_validation()
%MAIN_VALIDATION  Programme principal de validation numérique.
%
%  Ce script exécute les quatre cas de test décrits dans la Section 6
%  de l'article. Il peut être appelé directement :
%
%    >> main_validation

clear; clc; close all;

% Add the source tree to the MATLAB path (works from anywhere).
here = fileparts(mfilename('fullpath'));      % .../programmes/tests
root = fileparts(here);                        % .../programmes
addpath(genpath(fullfile(root,'src'))); addpath(here);

% Horodatage pour les logs
t_start_global = tic;
fprintf('===============================================================\n');
fprintf('  SEI ASYMPTOTIC VALIDATION PROGRAM\n');
fprintf('  Started : %s\n', char(datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss')));
fprintf('===============================================================\n');

% =====================================================================
%  Configuration des régimes à tester
% =====================================================================
% Pour les cas 1 et 2, on utilise le régime intermédiaire par défaut
% (le plus représentatif d'une SEI mûre)
regime_main = 'intermediate';

% =====================================================================
%  CAS 1 : Trilayer stationnaire
% =====================================================================
fprintf('\n>>> CONFIGURING PARAMETERS FOR REGIME : %s\n', regime_main);
P1 = config_parameters(regime_main);

results_case1 = test_case1_steady(P1);

% Génération de la figure de comparaison de champs
figures_dir = fullfile(pwd, 'figures');
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end
plot_field_comparison(results_case1.sol_full, results_case1.sol_eff, 'u', ...
    fullfile(figures_dir, sprintf('fig1_displacement_%s.png', P1.regime)));
plot_field_comparison(results_case1.sol_full, results_case1.sol_eff, 'c', ...
    fullfile(figures_dir, sprintf('fig1_concentration_%s.png', P1.regime)));
plot_field_comparison(results_case1.sol_full, results_case1.sol_eff, 'sigma', ...
    fullfile(figures_dir, sprintf('fig1_stress_%s.png', P1.regime)));

% =====================================================================
%  CAS 2 : Chargement cyclique transitoire
% =====================================================================
% Réduire la durée totale pour accélérer le test
P2 = config_parameters(regime_main);
P2.t_final = 0.5 * P2.T_cycle;  % un demi-cycle pour la démo
P2.dt      = 5.0;               % pas plus grand

results_case2 = test_case2_transient(P2);

% =====================================================================
%  CAS 3 : Endommagement et fissuration (régime critique)
% =====================================================================
P3 = config_parameters('critical');

results_case3 = test_case3_fracture(P3);

% Génération des courbes cohésives (Figure 3 de l'article)
plot_cohesive_curves(P3, fullfile(figures_dir, 'fig3_cohesive.png'));

% =====================================================================
%  CAS 4 : Étude de convergence sur 4 régimes (Figure 5)
% =====================================================================
conv_data = test_case4_convergence();

% =====================================================================
%  Bilan global
% =====================================================================
fprintf('\n===============================================================\n');
fprintf('  GLOBAL SUMMARY\n');
fprintf('===============================================================\n');

fprintf('\n  Case 1 (steady-state, %s)   : %s  (max err = %.2f%%)\n', ...
        P1.regime, ...
        ternary(results_case1.pass, 'PASS', 'FAIL'), ...
        100 * results_case1.err_max);

fprintf('  Case 2 (transient, %s)      : %s  (max err = %.2f%%)\n', ...
        P2.regime, ...
        ternary(results_case2.pass, 'PASS', 'FAIL'), ...
        100 * results_case2.err_max);

fprintf('  Case 3 (fracture, critical)  : %s  (max err = %.2f%%)\n', ...
        ternary(results_case3.pass, 'PASS', 'FAIL'), ...
        100 * results_case3.err_max);

fprintf('  Case 4 (convergence)         : %s  (O(eps) verified)\n', ...
        ternary(conv_data.pass, 'PASS', 'FAIL'));

fprintf('\n  Total elapsed time : %.1f s\n', toc(t_start_global));
fprintf('  End : %s\n', char(datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss')));
fprintf('===============================================================\n');

% Sauvegarde du bilan global
summary = struct(...
    'case1', results_case1, ...
    'case2', results_case2, ...
    'case3', results_case3, ...
    'case4', conv_data, ...
    'elapsed_time', toc(t_start_global));

results_dir = fullfile(pwd, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
save(fullfile(results_dir, 'summary_all_cases.mat'), 'summary', '-v7.3');
fprintf('\nGlobal summary saved to results/summary_all_cases.mat\n');

end  % <-- terminate main_validation (was missing: 'ternary' had an 'end' but
     %     the main function did not, which is an illegal mixed convention).

% =====================================================================
%  Fonction locale
% =====================================================================
function s = ternary(cond, val_true, val_false)
%TERNARY  Opérateur ternaire
if cond
    s = val_true;
else
    s = val_false;
end
end
