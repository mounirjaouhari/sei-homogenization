function plot_field_comparison(sol_full, sol_eff, field_name, output_path)
%PLOT_FIELD_COMPARISON  Compare visuellement les champs FULL et EFF.
%
%   plot_field_comparison(sol_full, sol_eff, field_name, output_path)
%   génère un graphique comparant le champ spécifié entre le modèle FULL
%   et le modèle EFFECTIF sur le domaine 1D.
%
%   Entrées :
%     sol_full    - solution du modèle FULL (solve_steady_state)
%     sol_eff     - solution du modèle EFF (solve_steady_state)
%     field_name  - nom du champ à comparer : 'u', 'c', ou 'sigma'
%     output_path - chemin de sauvegarde (optionnel)
%
%   Voir aussi : solve_steady_state, compute_l2_error

% Extraction des champs
field_full = sol_full.(field_name);
nodes_full = sol_full.nodes;
field_eff  = sol_eff.(field_name);
nodes_eff  = sol_eff.nodes;

% Conversion en μm pour l'affichage
x_full_um = nodes_full * 1e6;
x_eff_um  = nodes_eff  * 1e6;

% Unités selon le champ
switch field_name
    case 'u'
        y_label = 'Displacement u (nm)';
        y_full  = field_full * 1e9;  % m -> nm
        y_eff   = field_eff  * 1e9;
    case 'c'
        y_label = 'Concentration c (mol/m^3)';
        y_full  = field_full;
        y_eff   = field_eff;
    case 'sigma'
        y_label = 'Stress \sigma (MPa)';
        y_full  = field_full * 1e-6;  % Pa -> MPa
        y_eff   = field_eff  * 1e-6;
    otherwise
        error('plot_field_comparison:unknownField', ...
              'Champ inconnu : %s. Utiliser u/c/sigma.', field_name);
end

% Création de la figure
fig = figure('Position', [100, 100, 800, 500]);
plot(x_full_um, y_full, 'b-', 'LineWidth', 2, 'DisplayName', 'Full model');
hold on;
plot(x_eff_um, y_eff, 'r--', 'LineWidth', 2, 'DisplayName', 'Effective model');
hold off;

grid on;
xlabel('Position x (\mu m)', 'FontSize', 12);
ylabel(y_label, 'FontSize', 12);
title(sprintf('Comparison of %s field (%s regime)', ...
              field_name, sol_full.P.regime), 'FontSize', 13);
legend('Location', 'best', 'FontSize', 11);
set(gca, 'FontSize', 11);

% Ligne verticale à la position de l'interface (modèle FULL)
x_iface_um = sol_full.P.L * 1e6;
yl = ylim;                 % capturer le vecteur [ymin ymax] AVANT de l'utiliser
hold on;
plot([x_iface_um, x_iface_um], yl, 'k:', 'LineWidth', 1);
plot([x_iface_um + sol_full.P.eps * 1e6, x_iface_um + sol_full.P.eps * 1e6], ...
     yl, 'k:', 'LineWidth', 1);
text(x_iface_um + 0.5 * sol_full.P.eps * 1e6, yl(2) * 0.9, ...
     'SEI', 'HorizontalAlignment', 'center', 'FontSize', 10);
hold off;

% Sauvegarde
if nargin >= 4 && ~isempty(output_path)
    saveas(fig, output_path);
    fprintf('plot_field_comparison: figure saved to %s\n', output_path);
end
end
