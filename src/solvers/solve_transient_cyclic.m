function solution = solve_transient_cyclic(P, model_type, BC, options)
%SOLVE_TRANSIENT_CYCLIC  Résout le problème transitoire chémo-mécanique cyclique.
%
%   solution = solve_transient_cyclic(P, model_type, BC, options)
%   résout le problème couplé sur [0, t_final] avec un chargement
%   sinusoïdal représentant un cycle de charge/décharge (C/2).
%
%   Schéma temporel : Euler implicite
%     (M/dt + K) * c^{n+1} = M/dt * c^n + f_ext^{n+1}
%
%   Le champ de déplacement est résolu de manière quasi-statique à
%   chaque pas de temps (hypothèse de séparation des échelles).
%
%   Entrées :
%     P          - structure de paramètres
%     model_type - 'full' ou 'eff'
%     BC         - structure de CL (voir solve_steady_state)
%       .c_left_func  - fonction @(t) renvoyant c_left(t)
%       .c_right_func - fonction @(t) renvoyant c_right(t)
%       .u_left_func  - fonction @(t) renvoyant u_left(t)
%     options    - structure d'options :
%       .dt        - pas de temps (s) [défaut : P.dt]
%       .t_final   - durée totale (s) [défaut : P.t_final]
%       .verbose   - affichage de la progression (true/false)
%       .save_every- fréquence de sauvegarde (tous les N pas)
%
%   Sortie :
%     solution - structure :
%       .t        - (Nt x 1) vecteur temps
%       .u_hist   - (Nn x Nt) historique des déplacements
%       .c_hist   - (Nn x Nt) historique des concentrations
%       .sigma_hist - (Nn x Nt) historique des contraintes
%       .nodes, .elems, .regions, .model_type
%
%   Voir aussi : solve_steady_state

% =====================================================================
%  Options par défaut
% =====================================================================
if nargin < 4, options = struct(); end
if ~isfield(options, 'dt'),         options.dt = P.dt; end
if ~isfield(options, 't_final'),    options.t_final = P.t_final; end
if ~isfield(options, 'verbose'),    options.verbose = true; end
if ~isfield(options, 'save_every'), options.save_every = 1; end

% =====================================================================
%  Maillage et assemblage
% =====================================================================
switch lower(model_type)
    case 'full'
        [nodes, elems, regions] = build_trilayer_mesh(P);
        [K_u, ~, K_c, M_c] = assemble_full_model(nodes, elems, regions, P);
    case 'eff'
        [nodes, elems, regions] = build_eff_mesh(P);
        [K_u, ~, K_c, M_c] = assemble_eff_model(nodes, elems, regions, P);
    otherwise
        error('solve_transient_cyclic:unknownModel', ...
              'model_type doit être ''full'' ou ''eff''');
end

Nn = length(nodes);
dt = options.dt;
t_final = options.t_final;
t_vec = 0:dt:t_final;
Nt = length(t_vec);

% =====================================================================
%  Conditions initiales
% =====================================================================
% État initial : concentration nulle partout, sauf CL gauche
c = zeros(Nn, 1);
if isfield(BC, 'c_left_func')
    c(1) = BC.c_left_func(0);
end
if isfield(BC, 'c_right_func')
    c(Nn) = BC.c_right_func(0);
end

u = zeros(Nn, 1);
sigma = zeros(Nn, 1);

% Pré-allocation des historiques
u_hist     = zeros(Nn, Nt);
c_hist     = zeros(Nn, Nt);
sigma_hist = zeros(Nn, Nt);
u_hist(:, 1)     = u;
c_hist(:, 1)     = c;
sigma_hist(:, 1) = sigma;

% =====================================================================
%  Boucle temporelle (Euler implicite)
% =====================================================================
% Matrice constante pour la diffusion : A = M/dt + K
% (Les CL seront appliquées par pénalisation à chaque pas, modifiant A)
A_base = M_c / dt + K_c;

for n = 1:Nt-1
    t_np1 = t_vec(n+1);

    % ---------------------------------------------------------------
    %  1. Mise à jour de la concentration (diffusion implicite)
    % ---------------------------------------------------------------
    % Second membre : M/dt * c^n + f_ext^{n+1}
    b = (M_c / dt) * c;

    % Conditions aux limites dépendantes du temps
    c_dof = [];
    c_val = [];
    if isfield(BC, 'c_left_func')
        c_dof = [c_dof; 1];
        c_val = [c_val; BC.c_left_func(t_np1)];
    end
    if isfield(BC, 'c_right_func')
        c_dof = [c_dof; Nn];
        c_val = [c_val; BC.c_right_func(t_np1)];
    end

    % Application des CL (par pénalisation) - copie fraîche de A à chaque pas
    A_bc = A_base;
    b_bc = b;
    if ~isempty(c_dof)
        [A_bc, b_bc] = apply_dirichlet_bc(A_bc, b_bc, c_dof, c_val);
    end

    % Résolution
    c = A_bc \ b_bc;

    % ---------------------------------------------------------------
    %  2. Mise à jour du déplacement (quasi-statique)
    % ---------------------------------------------------------------
    f_u = compute_eigenstrain_force_transient(nodes, elems, regions, c, P, model_type);

    u_dof = [];
    u_val = [];
    if isfield(BC, 'u_left_func')
        u_dof = [u_dof; 1];
        u_val = [u_val; BC.u_left_func(t_np1)];
    end

    K_u_bc = K_u;
    f_u_bc = f_u;
    if ~isempty(u_dof)
        [K_u_bc, f_u_bc] = apply_dirichlet_bc(K_u_bc, f_u_bc, u_dof, u_val);
    end

    u = K_u_bc \ f_u_bc;

    % ---------------------------------------------------------------
    %  3. Calcul de la contrainte
    % ---------------------------------------------------------------
    sigma = compute_stress_transient(nodes, elems, regions, u, c, P);

    % ---------------------------------------------------------------
    %  Sauvegarde
    % ---------------------------------------------------------------
    u_hist(:, n+1)     = u;
    c_hist(:, n+1)     = c;
    sigma_hist(:, n+1) = sigma;

    % Affichage
    if options.verbose && mod(n, max(1, floor(Nt/10))) == 0
        fprintf('  t = %6.1f s / %6.1f s (%5.1f%%)\n', ...
                t_np1, t_final, 100 * n / (Nt - 1));
    end
end

% =====================================================================
%  Stockage de la solution
% =====================================================================
solution = struct(...
    't', t_vec', ...
    'u_hist', u_hist, ...
    'c_hist', c_hist, ...
    'sigma_hist', sigma_hist, ...
    'nodes', nodes, ...
    'elems', elems, ...
    'regions', regions, ...
    'model_type', lower(model_type), ...
    'P', P);
end

% ---------------------------------------------------------------------
function f_u = compute_eigenstrain_force_transient(nodes, elems, regions, c, P, model_type)
%COMPUTE_EIGENSTRAIN_FORCE_TRANSIENT  Force due à l'eigenstrain au temps courant

Nn = length(nodes);
f_u = zeros(Nn, 1);

for e = 1:size(elems, 1)
    idx = elems(e, :);
    x1 = nodes(idx(1));
    x2 = nodes(idx(2));

    region = regions.elems(e);
    switch region
        case 1
            C1d = P.C1d_1;  beta = P.beta1;
        case 2
            C1d = P.C1d_sei;  beta = P.beta_sei;
        case 3
            C1d = P.C1d_2;  beta = P.beta2;
    end

    c_elem = c(idx);
    c_avg  = 0.5 * (c_elem(1) + c_elem(2));
    fe = C1d * beta * (c_avg - P.c0) * [-1; 1];

    f_u(idx(1)) = f_u(idx(1)) + fe(1);
    f_u(idx(2)) = f_u(idx(2)) + fe(2);
end

if strcmp(model_type, 'eff')
    iL = regions.idx_interface_left;
    iR = regions.idx_interface_right;
    c_avg = 0.5 * (c(iL) + c(iR));
    f_int = P.K_eff * P.beta_eff * (c_avg - P.c0) * [-1; 1];
    f_u(iL) = f_u(iL) + f_int(1);
    f_u(iR) = f_u(iR) + f_int(2);
end
end

% ---------------------------------------------------------------------
function sigma = compute_stress_transient(nodes, elems, regions, u, c, P)
%COMPUTE_STRESS_TRANSIENT  Contrainte aux nœuds au temps courant

Nn = length(nodes);
sigma = zeros(Nn, 1);

for e = 1:size(elems, 1)
    idx = elems(e, :);
    x1 = nodes(idx(1));
    x2 = nodes(idx(2));
    L  = x2 - x1;

    region = regions.elems(e);
    switch region
        case 1
            C1d = P.C1d_1;  beta = P.beta1;
        case 2
            C1d = P.C1d_sei;  beta = P.beta_sei;
        case 3
            C1d = P.C1d_2;  beta = P.beta2;
    end

    du_dx = (u(idx(2)) - u(idx(1))) / L;
    c_avg = 0.5 * (c(idx(1)) + c(idx(2)));
    sigma_elem = C1d * (du_dx - beta * (c_avg - P.c0));

    sigma(idx(1)) = sigma(idx(1)) + 0.5 * sigma_elem;
    sigma(idx(2)) = sigma(idx(2)) + 0.5 * sigma_elem;
end
end
