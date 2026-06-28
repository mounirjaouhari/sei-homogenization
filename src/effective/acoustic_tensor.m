function Q = acoustic_tensor(C_tilde, n)
%ACOUSTIC_TENSOR  Tenseur acoustique pour matériau isotrope.
%
%   Q = acoustic_tensor(C_tilde, n) calcule le tenseur acoustique Q dans
%   la direction n pour un matériau élastique isotrope.
%
%   Le tenseur acoustique est défini par :
%     Q_ij = C_ijkl * n_k * n_l
%
%   Pour un matériau isotrope avec coefficients de Lamé (lambda, mu) :
%     Q = (lambda + mu) * n ⊗ n + mu * I
%
%   En 1D (n = [1]), Q est scalaire : Q = lambda + 2*mu
%
%   Entrées :
%     C_tilde - soit un scalaire (module 1D), soit une structure avec
%               champs .lambda et .mu (modules de Lamé 3D)
%     n       - vecteur unitaire (en 3D) ou scalaire 1 (en 1D)
%
%   Sortie :
%     Q - tenseur acoustique (scalaire en 1D, matrice 3x3 en 3D)
%
%   Voir aussi : effective_params

if isscalar(C_tilde)
    % Cas 1D : Q = C_tilde directement
    Q = C_tilde;
elseif isstruct(C_tilde) && isfield(C_tilde, 'lambda') && isfield(C_tilde, 'mu')
    % Cas 3D isotrope
    lambda = C_tilde.lambda;
    mu = C_tilde.mu;
    n = n(:)';
    Q = (lambda + mu) * (n' * n) + mu * eye(3);
else
    error('acoustic_tensor:invalidInput', ...
          'C_tilde doit être un scalaire (1D) ou une structure avec .lambda et .mu (3D)');
end
end
