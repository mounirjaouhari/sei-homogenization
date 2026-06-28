function [K_u, f_u, K_c, M_c] = assemble_eff_model(nodes, elems, regions, P)
%ASSEMBLE_EFF_MODEL  Assemblage des matrices du modèle EFFECTIF (interface).
%
%   [K_u, f_u, K_c, M_c] = assemble_eff_model(nodes, elems, regions, P)
%   assemble les matrices globales pour le modèle effectif où la SEI est
%   remplacée par une interface d'épaisseur nulle portant les conditions
%   de saut dérivées dans Theorem 3.1.
%
%   Conditions de saut (interface zéro-épaisseur) :
%     (1) Continuité de la traction : [[sigma * n]] = 0
%     (2) Loi d'interface chémo-mécanique :
%           sigma * n = K_eff * ([[u]] - beta_eff * (<c> - c0))
%     (3) Saut de concentration : [[c]] = R_eff * J
%     (4) Continuité du flux : [[J]] = 0,  J = -D * dc/dn
%
%   Implémentation : élément d'interface à 2 nœuds (nœud gauche électrode,
%   nœud droit électrolyte) avec matrices élémentaires :
%
%     Mécanique : K_int = K_eff * [1, -1; -1, 1]
%                 f_int = K_eff * beta_eff * (<c> - c0) * [-1; 1]
%
%     Diffusion : K_int = (1/R_eff) * [1, -1; -1, 1]
%                 M_int = 0 (pas de stockage dans l'interface)
%
%   Sorties :
%     K_u  - (Nn x Nn) matrice de raideur élastique globale
%     f_u  - (Nn x 1)  vecteur force
%     K_c  - (Nn x Nn) matrice de raideur diffusion globale
%     M_c  - (Nn x Nn) matrice de masse diffusion globale
%
%   Voir aussi : assemble_full_model, interface_element

Nn = length(nodes);
Ne = size(elems, 1);

% Initialisation
K_u = spalloc(Nn, Nn, 4 * Ne + 4);
K_c = spalloc(Nn, Nn, 4 * Ne + 4);
M_c = spalloc(Nn, Nn, 4 * Ne);
f_u = zeros(Nn, 1);

% Vecteur concentration initial
c_init = zeros(Nn, 1);

% =====================================================================
%  Boucle d'assemblage sur les éléments volumiques
% =====================================================================
for e = 1:Ne
    idx = elems(e, :);
    x1 = nodes(idx(1));
    x2 = nodes(idx(2));

    region = regions.elems(e);
    switch region
        case 1  % Électrode
            C1d  = P.C1d_1;
            beta = P.beta1;
            D    = P.D1;
        case 3  % Électrolyte
            C1d  = P.C1d_2;
            beta = P.beta2;
            D    = P.D2;
        otherwise
            error('assemble_eff_model:unknownRegion', ...
                  'Région inconnue : %d', region);
    end

    [Ke_u, fe_u] = elastic_element_stiff(x1, x2, C1d, beta, c_init(idx));
    [Ke_c, Me_c] = diffusion_element_stiff(x1, x2, D);

    for i = 1:2
        for j = 1:2
            K_u(idx(i), idx(j)) = K_u(idx(i), idx(j)) + Ke_u(i, j);
            K_c(idx(i), idx(j)) = K_c(idx(i), idx(j)) + Ke_c(i, j);
            M_c(idx(i), idx(j)) = M_c(idx(i), idx(j)) + Me_c(i, j);
        end
        f_u(idx(i)) = f_u(idx(i)) + fe_u(i);
    end
end

% =====================================================================
%  Ajout de l'élément d'interface (zéro épaisseur)
% =====================================================================
iL = regions.idx_interface_left;   % Nœud gauche de l'interface
iR = regions.idx_interface_right;  % Nœud droit de l'interface

% Élément d'interface mécanique
% K_int_u = K_eff * [1, -1; -1, 1]
K_int_u = P.K_eff * [1, -1; -1, 1];

% Force équivalente (couplage chémo-mécanique) :
% f_int = K_eff * beta_eff * (<c> - c0) * [-1; 1]
% avec <c> = (c(iL) + c(iR)) / 2
c_avg = 0.5 * (c_init(iL) + c_init(iR));
f_int_u = P.K_eff * P.beta_eff * (c_avg - 0) * [-1; 1];

% Assemblage
K_u(iL, iL) = K_u(iL, iL) + K_int_u(1, 1);
K_u(iL, iR) = K_u(iL, iR) + K_int_u(1, 2);
K_u(iR, iL) = K_u(iR, iL) + K_int_u(2, 1);
K_u(iR, iR) = K_u(iR, iR) + K_int_u(2, 2);

f_u(iL) = f_u(iL) + f_int_u(1);
f_u(iR) = f_u(iR) + f_int_u(2);

% Élément d'interface diffusion
% K_int_c = (1/R_eff) * [1, -1; -1, 1]
K_int_c = (1.0 / P.R_eff) * [1, -1; -1, 1];

K_c(iL, iL) = K_c(iL, iL) + K_int_c(1, 1);
K_c(iL, iR) = K_c(iL, iR) + K_int_c(1, 2);
K_c(iR, iL) = K_c(iR, iL) + K_int_c(2, 1);
K_c(iR, iR) = K_c(iR, iR) + K_int_c(2, 2);

% Note : pas de matrice de masse pour l'interface (capacité nulle)
end
