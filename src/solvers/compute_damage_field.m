function d = compute_damage_field(Delta, c_avg, P)
%COMPUTE_DAMAGE_FIELD  Champ d'endommagement pour la loi cohésive.
%
%   d = compute_damage_field(Delta, c_avg, P) retourne le champ
%   d'endommagement d ∈ [0, 1] associé à la séparation Delta et à la
%   concentration moyenne c_avg, en utilisant la loi cohésive
%   chémo-mécanique de Theorem 4.2.
%
%   Le calcul est délégué à chemo_cohesive_law, qui retourne également
%   la traction T et l'énergie Psi.
%
%   Entrées :
%     Delta  - (N x 1) séparations interfaciales (m)
%     c_avg  - (N x 1) ou scalaire, concentrations moyennes (mol/m^3)
%     P      - structure de paramètres
%
%   Sortie :
%     d - (N x 1) champ d'endommagement dans [0, 1]
%
%   Voir aussi : chemo_cohesive_law

[~, d, ~] = chemo_cohesive_law(Delta, c_avg, P);
end
