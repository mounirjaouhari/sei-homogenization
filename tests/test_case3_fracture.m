function results = test_case3_fracture(P)
%TEST_CASE3_FRACTURE  Cas de test 3 : endommagement et propagation de fissure.
%
%   results = test_case3_fracture(P) exécute le cas de test 3 de la
%   Section 6.4 de l'article : validation de la loi cohésive
%   chémo-mécanique dans le régime critique.
%
%   Comparaison entre :
%     - FULL : endommagement explicite dans la SEI maillée (gradient damage)
%     - EFF  : loi cohésive à l'interface zéro-épaisseur
%
%   Métriques :
%     - Énergie de fissuration prédite (formule) vs. numérique (intégrée)
%     - Longueur de fissure en fonction du chargement
%     - Profil de traction cohésive
%
%   Entrée :
%     P - structure de paramètres (régime critique recommandé)
%
%   Sortie :
%     results - structure avec champs :
%       .Gc_pred      - énergie prédite (formule eps * sum fi * gi)
%       .Gc_num_full  - énergie dissipée numérique (FULL)
%       .Gc_num_eff   - énergie dissipée numérique (EFF)
%       .err_Gc       - erreur relative sur Gc
%       .cohesive_data - courbes T(Delta) pour différents c_avg
%       .pass         - true si erreur < 10%
%
%   Voir aussi : chemo_cohesive_law, compute_damage_field

fprintf('\n===============================================================\n');
fprintf('  TEST CASE 3 : Fracture in critical regime\n');
fprintf('===============================================================\n');

% Forcer le régime critique si pas déjà
if ~strcmpi(P.regime, 'critical')
    warning('test_case3_fracture:notCritical', ...
            'Le régime actuel est %s ; le cas 3 est conçu pour le régime critique.', ...
            P.regime);
end

% =====================================================================
%  1. Calcul de l'énergie de fissuration prédite (Theorem 4.2)
% =====================================================================
Gc_pred = P.Gc;
fprintf('\n[1/4] Predicted fracture energy G_c = %.4e J/m^2\n', Gc_pred);
fprintf('      (eps_SEI = %.0f nm, <g> = %.2e J/m^3)\n', ...
        P.eps * 1e9, ...
        P.f_LiF * P.g_LiF + P.f_Li2CO3 * P.g_Li2CO3 + P.f_org * P.g_org);

% =====================================================================
%  2. Génération des courbes cohésives pour différents c_avg
% =====================================================================
fprintf('\n[2/4] Generating cohesive law curves...\n');
Delta_vec = linspace(0, 3 * P.delta_c, 200);
c_values  = [0, 0.1 * P.cmax, 0.2 * P.cmax];

cohesive_data = struct();
cohesive_data.Delta = Delta_vec;
cohesive_data.c_values = c_values;
cohesive_data.T = zeros(length(Delta_vec), length(c_values));
cohesive_data.d = zeros(length(Delta_vec), length(c_values));
cohesive_data.Psi = zeros(length(Delta_vec), length(c_values));

for k = 1:length(c_values)
    [T, d, Psi] = chemo_cohesive_law(Delta_vec, c_values(k), P);
    cohesive_data.T(:, k)   = T;
    cohesive_data.d(:, k)   = d;
    cohesive_data.Psi(:, k) = Psi;
end
fprintf('  Cohesive curves generated for %d concentrations\n', length(c_values));

% =====================================================================
%  3. Simulation quasi-statique avec propagation de fissure
%     Modèle EFF : interface cohésive, chargement monotone
% =====================================================================
fprintf('\n[3/4] Simulating crack propagation (effective model)...\n');
n_steps = 50;
Delta_applied = linspace(0, 2 * P.delta_c, n_steps);
T_response = zeros(n_steps, 1);
d_response = zeros(n_steps, 1);
energy_dissipated = zeros(n_steps, 1);
c_avg_test = 0;  % Pas de gonflement chimique pour ce test

for n = 1:n_steps
    [T, d, Psi] = chemo_cohesive_law(Delta_applied(n), c_avg_test, P);
    T_response(n) = T(1);
    d_response(n) = d(1);
    energy_dissipated(n) = Psi(1);
end

% Énergie dissipée totale = énergie à la séparation finale
% (intégrale de T d'Delta jusqu'à 2*delta_c)
Gc_num_eff = trapz(Delta_applied, T_response);
fprintf('  Numerical G_c (EFF) = %.4e J/m^2\n', Gc_num_eff);

% =====================================================================
%  4. Simulation FULL (gradient damage dans la SEI)
%     On résout un problème 1D dans la SEI avec endommagement
% =====================================================================
fprintf('\n[4/4] Simulating crack propagation (full model with damage)...\n');

% Maillage fin de la SEI
n_sei_fine = 50;
x_sei = linspace(0, P.eps, n_sei_fine + 1)';
L_elem = P.eps / n_sei_fine;

% Champ d'endommagement dans la SEI (initialisé à 0)
d_field = zeros(n_sei_fine + 1, 1);
u_field = zeros(n_sei_fine + 1, 1);

% Simulation : chargement monotone en déplacement aux bords
n_load_steps = 50;
Delta_load_full = linspace(0, 2 * P.delta_c, n_load_steps);
T_full = zeros(n_load_steps, 1);
d_full_avg = zeros(n_load_steps, 1);

for step = 1:n_load_steps
    % Déplacement imposé
    u_field(1) = 0;
    u_field(end) = Delta_load_full(step);

    % Itération de Newton pour résoudre le problème non-linéaire
    [u_field, d_field] = solve_damage_step(u_field, d_field, x_sei, P);

    % Traction = C1d * (1-d)^2 * du/dx à gauche
    du_dx_left = (u_field(2) - u_field(1)) / L_elem;
    d_field = d_field(:);                       % garantir un vecteur colonne
    T_full(step) = P.C1d_sei * (1 - d_field(1))^2 * du_dx_left;
    d_full_avg(step) = mean(d_field);
end

% Énergie dissipée (FULL) = intégrale de T d'Delta
Gc_num_full = trapz(Delta_load_full, T_full);
fprintf('  Numerical G_c (FULL) = %.4e J/m^2\n', Gc_num_full);

% =====================================================================
%  Calcul des erreurs
% =====================================================================
err_Gc_eff  = abs(Gc_pred - Gc_num_eff)  / Gc_pred;
err_Gc_full = abs(Gc_pred - Gc_num_full) / Gc_pred;

fprintf('\n---------------------------------------------------------------\n');
fprintf('  Fracture energy comparison :\n');
fprintf('    Predicted (formula)  :  G_c = %.4e J/m^2\n', Gc_pred);
fprintf('    Numerical (EFF)      :  G_c = %.4e J/m^2  (err = %.2f%%)\n', ...
        Gc_num_eff,  100 * err_Gc_eff);
fprintf('    Numerical (FULL)     :  G_c = %.4e J/m^2  (err = %.2f%%)\n', ...
        Gc_num_full, 100 * err_Gc_full);
fprintf('---------------------------------------------------------------\n');

err_max = max([err_Gc_eff, err_Gc_full]);
pass = err_max < P.error_accept;
if pass
    fprintf('  PASS : max error = %.2f%% < %.0f%% acceptance\n', ...
            100 * err_max, 100 * P.error_accept);
else
    fprintf('  FAIL : max error = %.2f%% >= %.0f%% acceptance\n', ...
            100 * err_max, 100 * P.error_accept);
end

% =====================================================================
%  Stockage
% =====================================================================
results = struct(...
    'Gc_pred',       Gc_pred, ...
    'Gc_num_full',   Gc_num_full, ...
    'Gc_num_eff',    Gc_num_eff, ...
    'err_Gc_full',   err_Gc_full, ...
    'err_Gc_eff',    err_Gc_eff, ...
    'err_max',       err_max, ...
    'pass',          pass, ...
    'cohesive_data', cohesive_data, ...
    'Delta_applied', Delta_applied, ...
    'T_response',    T_response, ...
    'd_response',    d_response, ...
    'Delta_load_full', Delta_load_full, ...
    'T_full',        T_full, ...
    'd_full_avg',    d_full_avg);

results_dir = fullfile(pwd, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
save(fullfile(results_dir, 'case3_fracture.mat'), 'results', '-v7.3');
fprintf('\nResults saved to results/case3_fracture.mat\n');
end

% ---------------------------------------------------------------------
function [u, d] = solve_damage_step(u, d, x, P)
%SOLVE_DAMAGE_STEP  Résout un pas non-linéaire du problème d'endommagement.
%
%   Modèle 1D simplifié dans la SEI :
%     sigma = C1d * (1-d)^2 * (du/dx - beta * c)
%     équilibrium : d(sigma)/dx = 0  =>  sigma = constant
%     d évolution : d = max(d_prev, max(0, (Y - Y0) / (Y1 - Y0)))
%     avec Y = sigma^2 / (2 * C1d * (1-d)^4) [taux de restitution]

n = length(x);
L_elem = x(2) - x(1);
c_avg = 0;  % pas de chimie ici

% Itérations de point fixe (simplifié)
for iter = 1:10
    % Calcul de sigma aux éléments (constant par hypothèse)
    sigma_const = P.C1d_sei * (1 - mean(d))^2 * ...
                  ((u(end) - u(1)) / P.eps - P.beta_sei * c_avg);

    % Champ de déplacement linéaire (équilibre => sigma constant)
    u = linspace(u(1), u(end), n)';

    % Mise à jour de l'endommagement
    Y = sigma_const^2 ./ (2 * P.C1d_sei * (1 - d + 1e-10).^4);  % ./ (élément par élément)
    Y0 = P.w1 / (2 * P.C1d_sei);
    d_new = max(d, max(0, min(1, (Y - Y0) / (P.w0 * P.w1))));

    % Test de convergence
    if norm(d_new - d) < 1e-8
        d = d_new;
        break;
    end
    d = d_new;
end
end
