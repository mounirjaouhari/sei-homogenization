function [K, f] = apply_dirichlet_bc(K, f, dof, value)
%APPLY_DIRICHLET_BC  Applique une condition de Dirichlet par pénalisation.
%
%   [K, f] = apply_dirichlet_bc(K, f, dof, value) impose u(dof) = value
%   en modifiant la matrice K et le vecteur f par la méthode de
%   pénalisation (grand nombre sur la diagonale).
%
%   Cette méthode préserve la structure creuse de K et évite le
%   réordonnancement des DDL. La valeur de pénalisation est choisie
%   grande devant les autres termes diagonaux.
%
%   Entrées :
%     K     - (N x N) matrice (creuse ou pleine)
%     f     - (N x 1) vecteur second membre
%     dof   - indice du DDL à contraindre (entier ou vecteur d'entiers)
%     value - valeur imposée (scalaire ou vecteur de même taille que dof)
%
%   Sorties :
%     K - matrice modifiée
%     f - vecteur modifié
%
%   Voir aussi : solve_steady_state, solve_transient_cyclic

N = length(f);
penalty = 1e20 * max(abs(diag(K)), [], 'all');
if penalty == 0
    penalty = 1e20;
end

% Conversion en vecteurs pour uniformiser
dof = dof(:);
value = value(:);
assert(length(dof) == length(value), ...
       'apply_dirichlet_bc:sizeMismatch', ...
       'dof et value doivent avoir la même taille');

for k = 1:length(dof)
    i = dof(k);
    if i < 1 || i > N
        error('apply_dirichlet_bc:outOfBounds', ...
              'DDL hors bornes : %d (N = %d)', i, N);
    end
    % Méthode de pénalisation
    f = f - K(:, i) * value(k);  % Report sur le second membre
    K(:, i) = 0;                  % Annulation de la colonne
    K(i, :) = 0;                  % Annulation de la ligne
    K(i, i) = penalty;            % Pénalisation diagonale
    f(i) = penalty * value(k);    % Second membre cohérent
end
end
