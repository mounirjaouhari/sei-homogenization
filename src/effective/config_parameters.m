function P = config_parameters(regime_name)
%CONFIG_PARAMETERS  Paramètres matériaux et exponent de scaling pour la SEI.
%
%   P = config_parameters(regime_name) retourne une structure P contenant
%   tous les paramètres physiques et numériques pour le régime asymptotique
%   spécifié.
%
%   Régimes supportés :
%     'soft'         - alpha = 1.5,  SEI très compliant (organique riche)
%     'intermediate' - alpha = 1.0,  SEI mixte standard
%     'critical'     - alpha = 0.0,  SEI inorganique, fissuration
%     'stiff'        - alpha = -0.5, SEI dense céramique
%
%   Paramètres retournés :
%     P.alpha, P.gamma, P.delta  - Exposants de scaling
%     P.eps   - Épaisseur SEI (m)
%     P.L     - Longueur macroscopique (m)
%     P.E1, P.E2, P.E_sei     - Modules de Young (Pa)
%     P.nu1, P.nu2, P.nu_sei  - Coefficients de Poisson
%     P.D1, P.D2, P.D_sei     - Diffusivités Li (m^2/s)
%     P.beta1, P.beta_sei     - Coefficients de Vegard (m^3/mol)
%     P.c0, P.cmax             - Concentrations de référence et max (mol/m^3)
%     P.K_eff                  - Raideur effective de l'interface (Pa/m)
%     P.beta_eff               - Couplage chémo-mécanique effectif (m^3/mol)
%     P.R_eff                  - Résistance de diffusion effective (s/m)
%     P.Gc                     - Énergie de fissuration (J/m^2)
%
%   Voir aussi : effective_params, build_trilayer_mesh

% =====================================================================
%  Sélection du régime
% =====================================================================
switch lower(regime_name)
    case 'soft'
        alpha = 1.5;  gamma = 1.5;  delta = 0.0;
    case 'intermediate'
        alpha = 1.0;  gamma = 1.2;  delta = 0.0;
    case 'critical'
        alpha = 0.0;  gamma = 1.0;  delta = 0.0;
    case 'stiff'
        alpha = -0.5; gamma = 0.8;  delta = 0.0;
    otherwise
        error('config_parameters:unknownRegime', ...
              'Régime inconnu : %s. Utiliser soft/intermediate/critical/stiff.', ...
              regime_name);
end

% =====================================================================
%  Géométrie
% =====================================================================
P.L        = 10e-6;          % Longueur caractéristique (10 μm)
P.eps      = 30e-9;          % Épaisseur SEI (30 nm par défaut)
P.h        = 1.0;            % Demi-épaisseur normalisée
P.regime   = lower(regime_name);
P.alpha    = alpha;
P.gamma    = gamma;
P.delta    = delta;

% =====================================================================
%  Propriétés électrode (graphite lithié)
% =====================================================================
P.E1       = 30e9;           % Module de Young graphite (Pa)
P.nu1      = 0.3;            % Poisson
P.D1       = 1e-12;          % Diffusivité Li dans graphite (m^2/s)
P.beta1    = 5e-3;           % Coefficient de Vegard graphite (m^3/mol)
P.rho1     = 2260;           % Masse volumique (kg/m^3)

% =====================================================================
%  Propriétés électrolyte (liquide, mécaniquement mou)
% =====================================================================
P.E2       = 1e6;            % Module très faible (Pa)
P.nu2      = 0.49;           % Proche incompressible
P.D2       = 1e-10;          % Diffusivité Li en électrolyte (m^2/s)
P.beta2    = 0.0;            % Pas de couplage chémo-mécanique

% =====================================================================
%  Propriétés SEI effectives (homogénéisées)
%  Valeurs typiques pour composition mixte (f_inorg ~ 0.6)
% =====================================================================
switch lower(regime_name)
    case 'soft'
        E_sei  = 2e9;     D_sei = 5e-15;   beta_sei = 1e-2;
    case 'intermediate'
        E_sei  = 15e9;    D_sei = 5e-15;   beta_sei = 5e-3;
    case 'critical'
        E_sei  = 25e9;    D_sei = 1e-15;   beta_sei = 3e-3;
    case 'stiff'
        E_sei  = 50e9;    D_sei = 5e-16;   beta_sei = 2e-3;
end
P.E_sei    = E_sei;
P.nu_sei   = 0.28;
P.D_sei    = D_sei;
P.beta_sei = beta_sei;

% Coefficient de Lamé pour le modèle 1D
P.lambda_sei = P.E_sei * P.nu_sei / ((1 + P.nu_sei) * (1 - 2 * P.nu_sei));
P.mu_sei     = P.E_sei / (2 * (1 + P.nu_sei));
P.C1d_sei    = P.lambda_sei + 2 * P.mu_sei;  % Module contrainte-déformation 1D

% De même pour l'électrode et l'électrolyte
P.lambda1 = P.E1 * P.nu1 / ((1 + P.nu1) * (1 - 2 * P.nu1));
P.mu1     = P.E1 / (2 * (1 + P.nu1));
P.C1d_1   = P.lambda1 + 2 * P.mu1;

P.lambda2 = P.E2 * P.nu2 / ((1 + P.nu2) * (1 - 2 * P.nu2));
P.mu2     = P.E2 / (2 * (1 + P.nu2));
P.C1d_2   = P.lambda2 + 2 * P.mu2;

% =====================================================================
%  Concentrations de référence
% =====================================================================
P.c0       = 0.0;            % Concentration de référence (mol/m^3)
P.cmax     = 3.0e4;          % Concentration max dans graphite (mol/m^3)

% =====================================================================
%  Paramètres d'endommagement (régime critique)
% =====================================================================
P.w0       = 1.0;            % Paramètre de régularisation d'endommagement
P.w1       = 1.0e6;          % Résistance à l'endommagement (J/m^3)
% NB : l'énergie volumique de fissuration g_i de chaque constituant est
% maintenant dérivée de sa résistance et de son module (g_i = sigma_r^2/(2E),
% J/m^3), voir compute_effective_params ci-dessous, plutôt que codée en dur.

% =====================================================================
%  Calcul des paramètres effectifs (Theorems 3.1 et 4.2)
% =====================================================================
P = compute_effective_params(P);

% =====================================================================
%  Paramètres temporels (pour cas transitoire)
% =====================================================================
P.T_cycle  = 3600;           % Période de cycle (1 h, C/2)
P.dt       = 1.0;            % Pas de temps (s)
P.t_final  = 2 * P.T_cycle;  % Durée totale : 2 cycles

% =====================================================================
%  Paramètres de maillage
% =====================================================================
P.n_elem_electrode   = 60;  % Éléments dans l'électrode
P.n_elem_electrolyte = 60;  % Éléments dans l'électrolyte
P.n_elem_sei         = 12;  % Éléments dans la SEI (modèle FULL)

% =====================================================================
%  Tolérances numériques
% =====================================================================
P.tol_linear    = 1e-10;    % Tolérance solveur linéaire
P.tol_newton    = 1e-8;     % Tolérance Newton (non-linéaire, cas 3)
P.max_newton    = 30;       % Itérations Newton max
P.error_accept  = 0.10;     % Critère d'acceptation (10%)

% Affichage
fprintf('config_parameters: régime = %s, alpha = %.2f, gamma = %.2f\n', ...
        P.regime, P.alpha, P.gamma);
fprintf('  K_eff = %.3e Pa/m, R_eff = %.3e s/m, G_c = %.3e J/m^2\n', ...
        P.K_eff, P.R_eff, P.Gc);
end

% ---------------------------------------------------------------------
function P = compute_effective_params(P)
%COMPUTE_EFFECTIVE_PARAMS  Calcule K_eff, beta_eff, R_eff, G_c (Theorems 3.1, 4.2)
%
%   À partir des exposants (alpha, gamma, delta) et des propriétés
%   homogénéisées de la SEI, calcule les paramètres effectifs de
%   l'interface limite.

% Coefficient de compliance 1D : A_tilde = 1 / C1d_sei
A_tilde = 1.0 / P.C1d_sei;

% Paramètres effectifs DIMENSIONNELLEMENT CORRECTS de l'interface, validés
% contre le modèle FULL (SEI explicitement maillée) par éléments finis :
% une couche de module C1d_sei, diffusivité D_sei et épaisseur t = h*eps agit
% comme un ressort de raideur C/t [Pa/m] et une résistance de diffusion t/D [s/m].
% (L'ancienne forme K_eff = (eps/L)^(alpha-1)*C1d/h était dimensionnellement Pa,
% pas Pa/m, et donnait une erreur L2 FULL/EFF de ~6e6 % ; voir matlab_val5.log.)
eta   = P.eps / P.L;                 % rapport d'échelle (info / régime)
t_sei = P.h * P.eps;                 % épaisseur (demi-)SEI dimensionnelle [m]

P.K_eff    = P.C1d_sei / t_sei;      % [Pa/m]  ressort à travers l'épaisseur
P.R_eff    = t_sei / P.D_sei;        % [s/m]   résistance de diffusion de la couche
P.beta_eff = P.beta_sei * t_sei;     % coefficient de saut de gonflement chimique
P.D_eff    = P.D_sei;                % diffusivité effective de l'interface
K_tilde    = P.K_eff;                % alias pour delta_c ci-dessous

% Énergie de fissuration effective (Théorème 4.2) :
%   G_c = eps_SEI * sum_i f_i * g_i,
% où g_i est une DENSITÉ VOLUMIQUE de fissuration [J/m^3] pour que
% [eps_SEI * g_i] = m * J/m^3 = J/m^2 (cohérence dimensionnelle).
% Choix documenté (reproductible) : g_i = sigma_{r,i}^2 / (2 E_i), la densité
% d'énergie élastique à la limite de résistance. (L'ancienne valeur g_i~1e9
% J/m^3 donnait G_c~85 J/m^2, ~280x la Table 3 de l'article : erreur d'unités.)
% Propriétés des constituants (résistance sigma_r [Pa], module E [Pa]) :
sig = struct('LiF', 250e6, 'Li2CO3', 200e6, 'org', 30e6);
Em  = struct('LiF', 65e9,  'Li2CO3', 60e9,  'org', 2e9);
P.f_LiF = 0.3;  P.f_Li2CO3 = 0.4;  P.f_org = 0.3;
P.g_LiF    = sig.LiF^2    / (2 * Em.LiF);     % J/m^3
P.g_Li2CO3 = sig.Li2CO3^2 / (2 * Em.Li2CO3);  % J/m^3
P.g_org    = sig.org^2    / (2 * Em.org);     % J/m^3
g_avg = P.f_LiF * P.g_LiF + P.f_Li2CO3 * P.g_Li2CO3 + P.f_org * P.g_org;
P.gc_volumetric = g_avg;     % densité volumique moyenne effective [J/m^3]
P.Gc = P.eps * g_avg;        % [J/m^2], dimensionnellement correct

% Séparation critique (Théorème 4.2) : delta_c = sqrt(2 G_c / K_tilde),
% avec K_tilde la raideur O(1) (et non K_eff), conformément à l'article.
P.delta_c = sqrt(2 * P.Gc / K_tilde);

end
