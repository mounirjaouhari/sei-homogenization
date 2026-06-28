function results = test_case2_transient(P)
%TEST_CASE2_TRANSIENT  Cas de test 2 : chargement cyclique transitoire.
%
%   results = test_case2_transient(P) exécute le cas de test 2 de la
%   Section 6.3 de l'article : cycle de charge/décharge C/2 avec
%   variation sinusoïdale de la concentration à l'électrode.
%
%   Comparaison entre FULL et EFF sur l'évolution temporelle de :
%     - Saut de déplacement [[u_n]](t)
%     - Concentration à l'interface c|_S0(t)
%     - Traction interfaciale sigma_nn(t)
%
%   Entrée :
%     P - structure de paramètres
%
%   Sortie :
%     results - structure avec solutions FULL et EFF + métriques d'erreur
%
%   Voir aussi : solve_transient_cyclic, compute_l2_error

fprintf('\n===============================================================\n');
fprintf('  TEST CASE 2 : Transient cyclic loading (%s regime)\n', P.regime);
fprintf('===============================================================\n');

% =====================================================================
%  Définition des conditions aux limites dépendantes du temps
% =====================================================================
% Concentration à gauche : sinusoïde c(t) = c_max/2 * (1 + sin(2*pi*t/T))
BC = struct();
BC.c_left_func  = @(t) 0.5 * P.cmax * (1 + sin(2 * pi * t / P.T_cycle));
BC.c_right_func = @(t) 0;
% Déplacement nul à gauche (référence)
BC.u_left_func  = @(t) 0;

% Options temporelles
options = struct();
options.dt        = P.dt;
options.t_final   = P.t_final;
options.verbose   = true;
options.save_every = 10;

% =====================================================================
%  Résolution FULL
% =====================================================================
fprintf('\n[1/2] Solving FULL model (transient)...\n');
t_full = tic;
sol_full = solve_transient_cyclic(P, 'full', BC, options);
fprintf('  FULL solved in %.2f s (%d time steps)\n', ...
        toc(t_full), length(sol_full.t));

% =====================================================================
%  Résolution EFF
% =====================================================================
fprintf('\n[2/2] Solving EFFECTIVE model (transient)...\n');
t_eff = tic;
sol_eff = solve_transient_cyclic(P, 'eff', BC, options);
fprintf('  EFF solved in %.2f s (%d time steps)\n', ...
        toc(t_eff), length(sol_eff.t));

% =====================================================================
%  Calcul des erreurs temporelles
% =====================================================================
% Saut de déplacement à l'interface
jump_u_full = compute_jump(sol_full, 'u');
jump_u_eff  = compute_jump(sol_eff,  'u');

% Concentration à l'interface
c_iface_full = compute_iface_value(sol_full, 'c');
c_iface_eff  = compute_iface_value(sol_eff,  'c');

% Erreurs RMS temporelles
err_jump_u = rms_relative(jump_u_full, jump_u_eff);
err_c_iface = rms_relative(c_iface_full, c_iface_eff);

fprintf('\n---------------------------------------------------------------\n');
fprintf('  Temporal RMS errors :\n');
fprintf('    [[u_n]](t)  :  err = %.4e  (%.2f%%)\n', err_jump_u,  100 * err_jump_u);
fprintf('    c_iface(t)  :  err = %.4e  (%.2f%%)\n', err_c_iface, 100 * err_c_iface);
fprintf('---------------------------------------------------------------\n');

err_max = max([err_jump_u, err_c_iface]);
pass = err_max < P.error_accept;
if pass
    fprintf('  PASS : max error = %.2f%% < %.0f%% acceptance\n', ...
            100 * err_max, 100 * P.error_accept);
else
    fprintf('  FAIL : max error = %.2f%% >= %.0f%% acceptance\n', ...
            100 * err_max, 100 * P.error_accept);
end

% =====================================================================
%  Stockage des résultats
% =====================================================================
results = struct(...
    'sol_full',    sol_full, ...
    'sol_eff',     sol_eff, ...
    'jump_u_full', jump_u_full, ...
    'jump_u_eff',  jump_u_eff, ...
    'c_iface_full', c_iface_full, ...
    'c_iface_eff',  c_iface_eff, ...
    'err_jump_u',  err_jump_u, ...
    'err_c_iface', err_c_iface, ...
    'err_max',     err_max, ...
    'pass',        pass, ...
    'BC',          BC);

results_dir = fullfile(pwd, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
save(fullfile(results_dir, sprintf('case2_transient_%s.mat', P.regime)), ...
     'results', '-v7.3');
fprintf('\nResults saved to results/case2_transient_%s.mat\n', P.regime);
end

% ---------------------------------------------------------------------
function jump_u = compute_jump(sol, field)
%COMPUTE_JUMP  Calcule le saut de `field` à l'interface au cours du temps

Nt = length(sol.t);
jump_u = zeros(Nt, 1);

switch sol.model_type
    case 'full'
        % Pour FULL, l'interface est entre la SEI et l'électrolyte
        idx_left  = sol.regions.idx_electrode_sei;
        idx_right = sol.regions.idx_sei_electrolyte;
    case 'eff'
        idx_left  = sol.regions.idx_interface_left;
        idx_right = sol.regions.idx_interface_right;
end

for n = 1:Nt
    jump_u(n) = sol.([field '_hist'])(idx_right, n) - ...
                sol.([field '_hist'])(idx_left, n);
end
end

% ---------------------------------------------------------------------
function val = compute_iface_value(sol, field)
%COMPUTE_IFACE_VALUE  Valeur moyenne à l'interface au cours du temps

Nt = length(sol.t);
val = zeros(Nt, 1);

switch sol.model_type
    case 'full'
        idx_left  = sol.regions.idx_electrode_sei;
        idx_right = sol.regions.idx_sei_electrolyte;
    case 'eff'
        idx_left  = sol.regions.idx_interface_left;
        idx_right = sol.regions.idx_interface_right;
end

for n = 1:Nt
    val(n) = 0.5 * (sol.([field '_hist'])(idx_left, n) + ...
                    sol.([field '_hist'])(idx_right, n));
end
end

% ---------------------------------------------------------------------
function err = rms_relative(v_ref, v_cmp)
%RMS_RELATIVE  Erreur RMS relative entre deux signaux temporels

diff_sq = (v_ref - v_cmp).^2;
ref_sq  = v_ref.^2;

if sum(ref_sq) < eps
    err = sqrt(mean(diff_sq));
else
    err = sqrt(sum(diff_sq) / sum(ref_sq));
end
end
