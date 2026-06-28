function plot_convergence_study(conv_data, output_path)
%PLOT_CONVERGENCE_STUDY  Trace les courbes de convergence O(eps).
%
%   plot_convergence_study(conv_data, output_path) génère la Figure 5
%   de l'article : convergence de l'erreur L2 relative entre le modèle
%   FULL et le modèle EFF en fonction de eps/L, pour les quatre régimes
%   asymptotiques.
%
%   Entrées :
%     conv_data   - structure avec champs :
%       .eps_over_L - (Neps x 1) valeurs de eps/L
%       .regimes    - cell array de noms de régime
%                     {'soft', 'intermediate', 'critical', 'stiff'}
%       .errors     - structure ou matrice (Neps x Nregimes) d'erreurs L2
%     output_path - chemin de sauvegarde de la figure (optionnel)
%
%   Voir aussi : test_case4_convergence, compute_l2_error

% Couleurs et marqueurs par régime
colors = {'blue', [0 0.6 0], 'red', [0.9 0.5 0]};
markers = {'o', 's', '^', 'd'};

% Création de la figure
fig = figure('Position', [100, 100, 800, 600]);
ax = axes('Parent', fig);
hold(ax, 'on');
grid(ax, 'on');

% Ligne de référence O(eps)
eps_vec = conv_data.eps_over_L;
ref_line = eps_vec / eps_vec(1) * conv_data.errors(1, 1) * 0.5;
loglog(ax, eps_vec, ref_line, 'k--', 'LineWidth', 1.5, 'DisplayName', 'O(\epsilon) reference');

% Tracé des courbes par régime
N_regimes = length(conv_data.regimes);
for r = 1:N_regimes
    err_r = conv_data.errors(:, r);
    loglog(ax, eps_vec, err_r, ...
           'Color', colors{r}, ...
           'Marker', markers{r}, ...
           'MarkerSize', 6, ...
           'LineWidth', 1.5, ...
           'DisplayName', sprintf('%s (\\alpha=%.1f)', ...
                                   conv_data.regimes{r}, ...
                                   get_alpha(conv_data.regimes{r})));
end

% Configuration des axes
xlabel(ax, 'Interphase thickness \epsilon/L (log scale)', 'FontSize', 12);
ylabel(ax, 'Relative error e_{rel} (log scale)', 'FontSize', 12);
title(ax, 'Convergence of the effective model vs. full FEM', 'FontSize', 13);
legend(ax, 'Location', 'NorthWest', 'FontSize', 10);
set(ax, 'FontSize', 11);
set(ax, 'XGrid', 'on', 'YGrid', 'on');

% Échelles logarithmiques (loglog) -- ne PAS forcer 'axis equal', qui repasse
% l'axe en échelle linéaire et écrase la pente de convergence.
set(ax, 'XScale', 'log', 'YScale', 'log');

% Sauvegarde
if nargin >= 2 && ~isempty(output_path)
    saveas(fig, output_path);
    fprintf('plot_convergence_study: figure saved to %s\n', output_path);
end
end

% ---------------------------------------------------------------------
function alpha = get_alpha(regime_name)
%GET_ALPHA  Exposant alpha associé au régime
switch lower(regime_name)
    case 'soft',         alpha = 1.5;
    case 'intermediate', alpha = 1.0;
    case 'critical',     alpha = 0.0;
    case 'stiff',        alpha = -0.5;
    otherwise,           alpha = NaN;
end
end
