function [K_int, f_int] = interface_element(P, c_left, c_right)
%INTERFACE_ELEMENT  Matrice élémentaire de l'interface chémo-mécanique.
%
%   [K_int, f_int] = interface_element(P, c_left, c_right) calcule la
%   matrice de raideur (mécanique) et le vecteur force de l'élément
%   d'interface zéro-épaisseur reliant le nœud gauche (côté électrode)
%   au nœud droit (côté électrolyte).
%
%   La loi d'interface (Theorem 3.1) est :
%
%     sigma * n = K_eff * ([[u]] - beta_eff * (<c> - c0))
%
%   où [[u]] = u_right - u_left, <c> = (c_left + c_right)/2.
%
%   La matrice de raideur élémentaire est :
%     K_int = K_eff * [1, -1; -1, 1]
%
%   Le vecteur force équivalent est :
%     f_int = K_eff * beta_eff * (<c> - c0) * [-1; 1]
%
%   Sorties :
%     K_int - (2x2) matrice de raideur de l'interface
%     f_int - (2x1) vecteur force équivalent
%
%   Voir aussi : assemble_eff_model, chemo_cohesive_law

% Validation des entrées
if P.K_eff <= 0
    error('interface_element:negativeStiffness', ...
          'K_eff doit être positif : K_eff = %e', P.K_eff);
end

% Matrice de raideur (mécanique)
K_int = P.K_eff * [1, -1; -1, 1];

% Force équivalente (eigenstrain chimique)
c_avg = 0.5 * (c_left + c_right);
f_int = P.K_eff * P.beta_eff * (c_avg - P.c0) * [-1; 1];
end
