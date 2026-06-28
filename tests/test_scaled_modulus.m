function out = test_scaled_modulus()
%TEST_SCALED_MODULUS  Convergence avec l'HYPOTHÈSE DE SCALING du module.
%
%   Teste la prédiction théorique centrale de l'article (régimes asymptotiques)
%   en appliquant au modèle FULL la loi d'échelle C^eps = (eps/L)^alpha * Ctilde,
%   c'est-à-dire un module de SEI qui DÉPEND de l'épaisseur (contrairement à
%   test_case4, qui utilise un module fixe). Le modèle EFF utilise le ressort
%   exact correspondant K_eff = C^eps/(h*eps) = (eps/L)^(alpha-1) * Ctilde/(h*L)... .
%
%   Question : le régime stiff (alpha<0) cesse-t-il de converger (O(1)) comme
%   l'affirme l'article, ou converge-t-il quand meme ? Réponse par simulation.
%
%   Sortie : out.eps_over_L, out.regimes, out.errors, out.rates

addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src')));

fprintf('\n===============================================================\n');
fprintf('  SCALED-MODULUS convergence  (C^eps = (eps/L)^alpha * Ctilde)\n');
fprintf('===============================================================\n');

eps_over_L = logspace(-4, log10(5e-2), 8);
Neps = length(eps_over_L);
regimes = {'soft', 'intermediate', 'critical', 'stiff'};
alphas  = [ 1.5,    1.0,           0.0,        -0.5];
Nreg = numel(regimes);
errors = zeros(Neps, Nreg);

eta_ref = 3e-3;   % référence où C^eps = module nominal du régime

for r = 1:Nreg
    P = config_parameters(regimes{r});
    alpha = alphas(r);
    C1d_nominal = P.C1d_sei;                      % module nominal du régime
    Ctilde = C1d_nominal / (eta_ref^alpha);       % => C^eps(eta_ref)=nominal
    fprintf('\n[%d/%d] %s (alpha=%.2f), Ctilde=%.3e Pa\n', r, Nreg, regimes{r}, alpha, Ctilde);

    for i = 1:Neps
        eta = eps_over_L(i);
        P.eps = eta * P.L;

        % MODULE SCALÉ du modèle FULL : C^eps = (eps/L)^alpha * Ctilde
        C_eps = (eta^alpha) * Ctilde;
        P.C1d_sei = C_eps;

        % EFF : ressort exact correspondant + résistance/gonflement de la couche
        t = P.h * P.eps;
        P.K_eff    = C_eps / t;          % = (eta^alpha) Ctilde / (h eps)
        P.R_eff    = t / P.D_sei;
        P.beta_eff = P.beta_sei * t;
        P.D_eff    = P.D_sei;

        BC = struct('u_left', 0, 'c_left', P.cmax, 'c_right', 0);
        try
            sf = solve_steady_state(P, 'full', BC);
            se = solve_steady_state(P, 'eff',  BC);
            errors(i, r) = compute_l2_error(sf.u, sf.nodes, se.u, se.nodes);
            fprintf('  eta=%.2e  C^eps=%.2e Pa  err=%.4e\n', eta, C_eps, errors(i, r));
        catch ME
            errors(i, r) = NaN;
            fprintf('  eta=%.2e  FAILED: %s\n', eta, ME.message);
        end
    end
end

fprintf('\n---------------------------------------------------------------\n');
fprintf('  Convergence rates (scaled modulus) :\n');
rates = nan(1, Nreg);
for r = 1:Nreg
    v = ~isnan(errors(:, r)) & errors(:, r) > 0;
    if sum(v) >= 2
        p = polyfit(log(eps_over_L(v)), log(errors(v, r)'), 1);
        rates(r) = p(1);
        fprintf('  %-15s (alpha=%+.1f) : rate = %.2f\n', regimes{r}, alphas(r), rates(r));
    else
        fprintf('  %-15s : insufficient/failed data\n', regimes{r});
    end
end
fprintf('---------------------------------------------------------------\n');

out = struct('eps_over_L', eps_over_L, 'regimes', {regimes}, ...
             'alphas', alphas, 'errors', errors, 'rates', rates);
end
