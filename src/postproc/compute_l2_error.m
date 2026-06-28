function err = compute_l2_error(field_full, nodes_full, field_eff, nodes_eff)
%COMPUTE_L2_ERROR  Erreur L2 relative entre deux champs 1D.
%
%   err = compute_l2_error(field_full, nodes_full, field_eff, nodes_eff)
%   calcule l'erreur L2 relative entre un champ de référence (modèle FULL)
%   et un champ à comparer (modèle EFF), sur leur maillage respectif.
%
%   L'erreur est définie comme :
%     err = ||field_full - field_eff_interp||_L2 / ||field_full||_L2
%
%   où field_eff_interp est l'interpolation du champ EFF sur le maillage FULL.
%
%   Entrées :
%     field_full - (N_full x 1) champ de référence
%     nodes_full - (N_full x 1) coordonnées du maillage FULL
%     field_eff  - (N_eff x 1) champ à comparer
%     nodes_eff  - (N_eff x 1) coordonnées du maillage EFF
%
%   Sortie :
%     err - erreur L2 relative (scalaire)
%
%   Voir aussi : solve_steady_state, plot_field_comparison

% Vérification des entrées
assert(length(field_full) == length(nodes_full), ...
       'compute_l2_error:sizeMismatch', 'field_full et nodes_full doivent avoir la même taille');
assert(length(field_eff) == length(nodes_eff), ...
       'compute_l2_error:sizeMismatch', 'field_eff et nodes_eff doivent avoir la même taille');

% Interpolation du champ EFF sur le maillage FULL.
% interp1 exige des abscisses STRICTEMENT croissantes ; or le maillage effectif
% comporte un nœud dédoublé à l'interface (épaisseur nulle) -> deux abscisses
% identiques. On trie puis on décale infinitésimalement les doublons pour
% conserver les DEUX valeurs (le saut) tout en rendant l'abscisse monotone.
x_src = nodes_eff(:);
f_src = field_eff(:);
[x_src, isrt] = sort(x_src);
f_src = f_src(isrt);
u_unique = unique(x_src);
if numel(u_unique) > 1
    tiny = 1e-6 * min(diff(u_unique));
else
    tiny = 1e-12;
end
for k = 2:numel(x_src)
    if x_src(k) <= x_src(k-1)
        x_src(k) = x_src(k-1) + tiny;
    end
end
field_eff_interp = interp1(x_src, f_src, nodes_full, 'linear', 'extrap');

% Calcul de la norme L2 par intégration trapézoïdale
% ||f||_L2 = sqrt(int f^2 dx)
norm_full_sq = trapz(nodes_full, field_full.^2);
norm_diff_sq = trapz(nodes_full, (field_full - field_eff_interp).^2);

if norm_full_sq < eps
    warning('compute_l2_error:zeroNorm', ...
            'Norme L2 du champ FULL nulle, retour de la norme absolue.');
    err = sqrt(norm_diff_sq);
else
    err = sqrt(norm_diff_sq / norm_full_sq);
end
end
