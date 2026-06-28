function [K_u, f_u, K_c, M_c] = assemble_full_model(nodes, elems, regions, P)
%ASSEMBLE_FULL_MODEL  Assemblage des matrices du modèle FULL (SEI maillée).
%
%   [K_u, f_u, K_c, M_c] = assemble_full_model(nodes, elems, regions, P)
%   assemble les matrices globales pour le problème chémo-mécanique 1D
%   avec la SEI explicitement maillée.
%
%   Le problème est découplé à l'itération (résolution séquentielle) :
%     1. Mécanique : K_u * u = f_u(u_mech, c)
%     2. Diffusion : (M_c/dt + K_c) * c^{n+1} = M_c/dt * c^n + f_ext
%
%   Sorties :
%     K_u  - (Nn x Nn) matrice de raideur élastique globale
%     f_u  - (Nn x 1)  vecteur force (eigenstrain, calculé pour c courant)
%     K_c  - (Nn x Nn) matrice de raideur diffusion globale
%     M_c  - (Nn x Nn) matrice de masse diffusion globale
%
%   Voir aussi : assemble_eff_model, elastic_element_stiff, diffusion_element_stiff

Nn = length(nodes);
Ne = size(elems, 1);

% Initialisation des matrices creuses
K_u = spalloc(Nn, Nn, 4 * Ne);
K_c = spalloc(Nn, Nn, 4 * Ne);
M_c = spalloc(Nn, Nn, 4 * Ne);
f_u = zeros(Nn, 1);

% Vecteur concentration (initial, pour le calcul de l'eigenstrain)
c_init = zeros(Nn, 1);

% Boucle d'assemblage
for e = 1:Ne
    % Nœuds de l'élément
    idx = elems(e, :);
    x1 = nodes(idx(1));
    x2 = nodes(idx(2));
    L  = x2 - x1;

    % Sélection des propriétés selon la région
    region = regions.elems(e);
    switch region
        case 1  % Électrode
            C1d  = P.C1d_1;
            beta = P.beta1;
            D    = P.D1;
        case 2  % SEI
            % La SEI a un module effectif dépendant de eps (scaling alpha)
            % C_sei = eps^alpha * C_tilde, où C_tilde est O(1)
            % Ici, P.E_sei est déjà la valeur effective pour le eps courant
            C1d  = P.C1d_sei;
            beta = P.beta_sei;
            D    = P.D_sei;
        case 3  % Électrolyte
            C1d  = P.C1d_2;
            beta = P.beta2;
            D    = P.D2;
        otherwise
            error('assemble_full_model:unknownRegion', ...
                  'Région inconnue : %d', region);
    end

    % Matrices élémentaires
    [Ke_u, fe_u] = elastic_element_stiff(x1, x2, C1d, beta, c_init(idx));
    [Ke_c, Me_c] = diffusion_element_stiff(x1, x2, D);

    % Assemblage
    for i = 1:2
        for j = 1:2
            K_u(idx(i), idx(j)) = K_u(idx(i), idx(j)) + Ke_u(i, j);
            K_c(idx(i), idx(j)) = K_c(idx(i), idx(j)) + Ke_c(i, j);
            M_c(idx(i), idx(j)) = M_c(idx(i), idx(j)) + Me_c(i, j);
        end
        f_u(idx(i)) = f_u(idx(i)) + fe_u(i);
    end
end

% Vérification de la symétrie
assert(norm(K_u - K_u', 'fro') / norm(K_u, 'fro') < 1e-12, ...
       'assemble_full_model:asymmetry', 'K_u n''est pas symétrique');
end
