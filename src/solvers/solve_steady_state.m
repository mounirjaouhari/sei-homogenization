function solution = solve_steady_state(P, model_type, BC)
%SOLVE_STEADY_STATE  Résout le problème stationnaire chémo-mécanique.
%
%   solution = solve_steady_state(P, model_type, BC) résout le problème
%   stationnaire couplé élastique-diffusion pour le modèle spécifié.
%
%   Le problème stationnaire se découple en deux sous-problèmes :
%     1. Diffusion : K_c * c = f_c  (régime permanent : dc/dt = 0)
%     2. Mécanique : K_u * u = f_u(c)  (eigenstrain calculé avec c résolu)
%
%   Entrées :
%     P          - structure de paramètres (config_parameters)
%     model_type - 'full' ou 'eff'
%     BC         - structure de conditions aux limites :
%       .u_left  - déplacement imposé à gauche (m)
%       .u_right - déplacement imposé à droite (m) [optionnel]
%       .c_left  - concentration imposée à gauche (mol/m^3)
%       .c_right - concentration imposée à droite (mol/m^3)
%       .J_right - flux imposé à droite (mol/(m^2.s)) [optionnel]
%
%   Sortie :
%     solution - structure :
%       .u       - (Nn x 1) champ de déplacement
%       .c       - (Nn x 1) champ de concentration
%       .sigma   - (Nn x 1) champ de contrainte
%       .nodes   - (Nn x 1) coordonnées
%       .regions - structure des régions
%       .model_type
%
%   Voir aussi : solve_transient_cyclic, assemble_full_model, assemble_eff_model

% =====================================================================
%  Construction du maillage
% =====================================================================
switch lower(model_type)
    case 'full'
        [nodes, elems, regions] = build_trilayer_mesh(P);
    case 'eff'
        [nodes, elems, regions] = build_eff_mesh(P);
    otherwise
        error('solve_steady_state:unknownModel', ...
              'model_type doit être ''full'' ou ''eff''');
end

Nn = length(nodes);

% =====================================================================
%  Assemblage des matrices
% =====================================================================
switch lower(model_type)
    case 'full'
        [K_u, f_u_init, K_c, M_c] = assemble_full_model(nodes, elems, regions, P);
    case 'eff'
        [K_u, f_u_init, K_c, M_c] = assemble_eff_model(nodes, elems, regions, P);
end

% =====================================================================
%  1. Résolution du problème de DIFFUSION stationnaire
%     K_c * c = f_c
%  =====================================================================
f_c = zeros(Nn, 1);

% Conditions aux limites pour la concentration
c_dof = [];
c_val = [];

if isfield(BC, 'c_left') && ~isempty(BC.c_left)
    c_dof = [c_dof; 1];
    c_val = [c_val; BC.c_left];
end
if isfield(BC, 'c_right') && ~isempty(BC.c_right)
    c_dof = [c_dof; Nn];
    c_val = [c_val; BC.c_right];
end

% Application des CL de Dirichlet
K_c_bc = K_c;
f_c_bc = f_c;
if ~isempty(c_dof)
    [K_c_bc, f_c_bc] = apply_dirichlet_bc(K_c_bc, f_c_bc, c_dof, c_val);
end

% Flux imposé à droite (Neumann) : ajout direct au second membre
if isfield(BC, 'J_right') && ~isempty(BC.J_right)
    f_c_bc(Nn) = f_c_bc(Nn) + BC.J_right;
end

% Résolution
c = K_c_bc \ f_c_bc;

% =====================================================================
%  2. Résolution du problème MÉCANIQUE avec eigenstrain calculé
%     K_u * u = f_u(c)
%  =====================================================================
% Recalcul de f_u avec le champ c résolu
f_u = compute_eigenstrain_force(nodes, elems, regions, c, P, model_type);

% Conditions aux limites pour le déplacement
u_dof = [];
u_val = [];

if isfield(BC, 'u_left') && ~isempty(BC.u_left)
    u_dof = [u_dof; 1];
    u_val = [u_val; BC.u_left];
end
if isfield(BC, 'u_right') && ~isempty(BC.u_right)
    u_dof = [u_dof; Nn];
    u_val = [u_val; BC.u_right];
end

K_u_bc = K_u;
f_u_bc = f_u;
if ~isempty(u_dof)
    [K_u_bc, f_u_bc] = apply_dirichlet_bc(K_u_bc, f_u_bc, u_dof, u_val);
end

% Résolution
u = K_u_bc \ f_u_bc;

% =====================================================================
%  Post-traitement : contrainte
% =====================================================================
sigma = compute_stress(nodes, elems, regions, u, c, P);

% =====================================================================
%  Stockage de la solution
% =====================================================================
solution = struct(...
    'u', u, ...
    'c', c, ...
    'sigma', sigma, ...
    'nodes', nodes, ...
    'elems', elems, ...
    'regions', regions, ...
    'model_type', lower(model_type), ...
    'P', P);
end

% ---------------------------------------------------------------------
function f_u = compute_eigenstrain_force(nodes, elems, regions, c, P, model_type)
%COMPUTE_EIGENSTRAIN_FORCE  Recalcule le vecteur force avec le champ c courant.
%
%   Pour chaque élément, la contribution de l'eigenstrain est :
%     fe = C1d * beta * (c_avg - c0) * [-1; 1]
%   avec c_avg = (c_i1 + c_i2)/2

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

% Contribution de l'interface (modèle effectif)
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
function sigma = compute_stress(nodes, elems, regions, u, c, P)
%COMPUTE_STRESS  Calcule le champ de contrainte aux nœuds.
%
%   sigma = C1d * (du/dx - beta * (c - c0))
%   Pour un élément linéaire : du/dx = (u2 - u1)/L (constant dans l'élément)

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

    % Attribution aux nœuds (moyenne pondérée)
    sigma(idx(1)) = sigma(idx(1)) + 0.5 * sigma_elem;
    sigma(idx(2)) = sigma(idx(2)) + 0.5 * sigma_elem;
end
end
