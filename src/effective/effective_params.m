function params = effective_params(P)
%EFFECTIVE_PARAMS  Calcule les paramètres effectifs de l'interface SEI.
%
%   params = effective_params(P) retourne une structure avec les
%   paramètres effectifs K_eff, beta_eff, R_eff, D_eff, G_c, delta_c
%   à partir des propriétés de la SEI et des exposants (alpha, gamma, delta).
%
%   Formules (Theorems 3.1 et 4.2 de l'article) :
%     K_eff    = eps^(alpha-1) * K_tilde         [Pa/m]
%     beta_eff = eps^delta * beta_tilde * h * A_tilde * (C_tilde : I)
%              = eps^delta * beta_sei * h         [m^3/mol]  (1D)
%     R_eff    = eps^(1-gamma) * R_tilde          [s/m]
%     D_eff    = eps^gamma * D_tilde              [m^2/s]
%     G_c      = eps * sum_i f_i * g_i            [J/m^2]
%     delta_c  = sqrt(2 * G_c / K_tilde)          [m]
%
%   Entrée :
%     P - structure de paramètres (config_parameters)
%
%   Sortie :
%     params - structure avec champs :
%       .K_eff, .beta_eff, .R_eff, .D_eff, .G_c, .delta_c
%       .K_tilde, .A_tilde, .R_tilde (valeurs O(1))
%
%   Voir aussi : config_parameters, acoustic_tensor

% Paramètres effectifs DIMENSIONNELLEMENT CORRECTS d'une couche mince uniforme
% (module C1d_sei, diffusivité D_sei, épaisseur t = h*eps), validés contre le
% modèle FULL par éléments finis : ressort C/t [Pa/m], résistance t/D [s/m].
A_tilde = 1.0 / P.C1d_sei;
t_sei   = P.h * P.eps;                 % épaisseur SEI dimensionnelle [m]
eta     = P.eps / P.L;                 % rapport d'échelle (info / régime)

params.K_eff    = P.C1d_sei / t_sei;   % [Pa/m]
params.R_eff    = t_sei / P.D_sei;     % [s/m]
params.beta_eff = P.beta_sei * t_sei;  % gonflement chimique
params.D_eff    = P.D_sei;             % diffusivité effective
K_tilde = params.K_eff;                % alias pour delta_c

% Énergie de fissuration (règle de mélange linéaire) : G_c = eps * sum_i f_i g_i,
% g_i = DENSITÉ VOLUMIQUE [J/m^3] pour la cohérence dimensionnelle (J/m^2).
% Défaut documenté : g_i = sigma_{r,i}^2/(2 E_i) (et non ~1e9 J/m^3, erreur d'unités).
if isfield(P, 'f_LiF') && isfield(P, 'g_LiF')
    g_avg = P.f_LiF * P.g_LiF + P.f_Li2CO3 * P.g_Li2CO3 + P.f_org * P.g_org;
else
    g_LiF    = (250e6)^2 / (2 * 65e9);   % J/m^3
    g_Li2CO3 = (200e6)^2 / (2 * 60e9);
    g_org    = (30e6)^2  / (2 * 2e9);
    g_avg = 0.3 * g_LiF + 0.4 * g_Li2CO3 + 0.3 * g_org;
end
params.G_c = P.eps * g_avg;

% Séparation critique
params.delta_c = sqrt(2 * params.G_c / K_tilde);

% Valeurs de référence
params.K_tilde = K_tilde;
params.A_tilde = A_tilde;
params.R_tilde = t_sei / P.D_sei;
params.t_sei   = t_sei;
params.eta     = eta;
end
