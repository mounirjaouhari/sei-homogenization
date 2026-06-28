function [T, d, Psi] = chemo_cohesive_law(Delta, c_avg, P)
%CHEMO_COHESIVE_LAW  Loi cohésive chémo-mécanique du régime critique (Thm 4.2).
%
%   [T, d, Psi] = chemo_cohesive_law(Delta, c_avg, P) renvoie, pour un vecteur
%   de séparations Delta [m], une concentration moyenne c_avg [mol/m^3] et la
%   structure de paramètres P :
%
%     T   - traction cohésive [Pa]      : T = (1-d)^2 * K_eff * Delta_eff
%     d   - endommagement [0,1]
%     Psi - énergie stockée+dissipée par unité de surface [J/m^2]
%
%   Séparation effective chémo-mécanique (le gonflement pré-déforme l'interface) :
%     Delta_eff = Delta - beta_eff * (c_avg - c0)
%
%   La loi est à adoucissement linéaire, calibrée pour que l'aire sous T(Delta)
%   vaille exactement l'énergie de fissuration G_c = P.Gc et que T s'annule à la
%   séparation critique delta_c = sqrt(2 G_c / K_eff) :
%     T(Delta_eff) = T_peak * (1 - Delta_eff/delta_c),  T_peak = sqrt(2 G_c K_eff)
%
%   Voir aussi : effective_params, test_case3_fracture

K  = P.K_eff;                 % raideur effective normale [Pa/m]
Gc = P.Gc;                    % énergie de fissuration [J/m^2]
c0 = P.c0;

if K <= 0 || Gc <= 0
    error('chemo_cohesive_law:badParams', 'K_eff et Gc doivent être > 0.');
end

delta_c = sqrt(2 * Gc / K);   % séparation critique [m]
T_peak  = sqrt(2 * Gc * K);   % = K * delta_c [Pa]

Delta = Delta(:);
Delta_eff = Delta - P.beta_eff * (c_avg - c0);

T   = zeros(size(Delta));
d   = zeros(size(Delta));
Psi = zeros(size(Delta));

for i = 1:numel(Delta)
    de = Delta_eff(i);
    if de <= 0
        % Pas de chargement effectif : interface intacte, traction nulle.
        T(i) = 0; d(i) = 0; Psi(i) = 0;
    elseif de >= delta_c
        % Au-delà de la séparation critique : rompu.
        T(i) = 0; d(i) = 1; Psi(i) = Gc;
    else
        % Adoucissement linéaire.
        T(i) = T_peak * (1 - de / delta_c);
        % Endommagement déduit de T = (1-d)^2 * K * de.
        ratio = T(i) / (K * de);
        ratio = min(max(ratio, 0), 1);
        d(i) = 1 - sqrt(ratio);
        % Aire sous la courbe jusqu'à de (énergie stockée + dissipée).
        Psi(i) = T_peak * de * (1 - de / (2 * delta_c));
    end
end
end
