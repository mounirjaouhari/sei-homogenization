function plot_cohesive_curves(P, output_path)
%PLOT_COHESIVE_CURVES  Trace les lois cohesives chémo-mécaniques.
%
%   plot_cohesive_curves(P, output_path) génère la Figure 3 de l'article :
%   la traction cohésive T en fonction de la séparation effective
%   Delta_eff, pour trois concentrations moyennes différentes.
%
%   Entrées :
%     P           - structure de paramètres (régime critique)
%     output_path - chemin de sauvegarde (optionnel)
%
%   Voir aussi : chemo_cohesive_law, compute_damage_field

% Plage de séparations (en nm)
Delta_nm = linspace(0, 5, 200);
Delta    = Delta_nm * 1e-9;  % m

% Trois concentrations moyennes
c_values = [0, 0.1 * P.cmax, 0.2 * P.cmax];
colors   = {'blue', 'red', [0.9 0.5 0]};
labels   = {'c = c_0 (no swelling)', ...
             'c = c_0 + 0.1 c_{max}', ...
             'c = c_0 + 0.2 c_{max}'};

% Création de la figure
fig = figure('Position', [100, 100, 800, 600]);
hold on;
grid on;

for k = 1:length(c_values)
    c_avg = c_values(k);
    [T, ~, ~] = chemo_cohesive_law(Delta, c_avg, P);
    T_MPa = T * 1e-6;  % Pa -> MPa
    plot(Delta_nm, T_MPa, 'Color', colors{k}, 'LineWidth', 2, ...
         'DisplayName', labels{k});
end

xlabel('Effective separation \Delta_{eff} (nm)', 'FontSize', 12);
ylabel('Cohesive traction T (MPa)', 'FontSize', 12);
title('Chemo-mechanical cohesive law at varying \langle c \rangle', 'FontSize', 13);
legend('Location', 'NorthEast', 'FontSize', 11);
set(gca, 'FontSize', 11);

% Marquer la séparation critique
delta_c_nm = P.delta_c * 1e9;
hold on;
plot([delta_c_nm, delta_c_nm], ylim, 'k:', 'LineWidth', 1);
text(delta_c_nm, max(ylim) * 0.9, ...
     sprintf('\\delta_c = %.2f nm', delta_c_nm), ...
     'HorizontalAlignment', 'right', 'FontSize', 10);
hold off;

% Sauvegarde
if nargin >= 2 && ~isempty(output_path)
    saveas(fig, output_path);
    fprintf('plot_cohesive_curves: figure saved to %s\n', output_path);
end
end
