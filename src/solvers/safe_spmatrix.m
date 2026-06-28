function S = safe_spmatrix(rows, cols, vals, N, M)
%SAFE_SPMATRIX  Construction robuste d'une matrice creuse.
%
%   S = safe_spmatrix(rows, cols, vals, N, M) construit une matrice
%   creuse N x M à partir des triplets (rows, cols, vals) en gérant les
%   doublons par sommation.
%
%   Cette fonction est un wrapper autour de sparse() qui assure que les
%   entrées dupliquées sont sommées (comportement natif de sparse), avec
%   vérifications de cohérence.
%
%   Entrées :
%     rows - (K x 1) indices de ligne
%     cols - (K x 1) indices de colonne
%     vals - (K x 1) valeurs
%     N    - nombre de lignes
%     M    - nombre de colonnes
%
%   Sortie :
%     S - (N x M) matrice creuse
%
%   Voir aussi : sparse

% Validation des entrées
K = length(rows);
assert(length(cols) == K && length(vals) == K, ...
       'safe_spmatrix:sizeMismatch', ...
       'rows, cols, vals doivent avoir la même taille');

assert(all(rows >= 1) && all(rows <= N) && all(cols >= 1) && all(cols <= M), ...
       'safe_spmatrix:outOfBounds', ...
       'Indices hors bornes');

% Construction de la matrice creuse
S = sparse(rows, cols, vals, N, M);
end
