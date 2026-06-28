function [Ke, Me] = diffusion_element_stiff(x1, x2, D)
%DIFFUSION_ELEMENT_STIFF  Matrices de raideur et de masse pour diffusion 1D.
%
%   [Ke, Me] = diffusion_element_stiff(x1, x2, D) calcule :
%     - Ke : matrice de raideur diffusion (2x2)
%            K_ij = int_Omega D * dN_i/dx * dN_j/dx dV
%     - Me : matrice de masse (2x2)
%            M_ij = int_Omega N_i * N_j dV
%
%   L'équation de diffusion transitoire est :
%     dc/dt = div(D grad c)
%   Formulation faible :
%     int_Omega dc/dt * v dV + int_Omega D * grad c * grad v dV = int_Gamma J*v dS
%
%   Discrétisation (Euler implicite) :
%     (M/dt + K) * c^{n+1} = M/dt * c^n + f_ext
%
%   Pour un élément linéaire 1D à 2 nœuds :
%     Ke = (D/L) * [[1, -1]; [-1, 1]]
%     Me = (L/6) * [[2, 1]; [1, 2]]
%
%   Sorties :
%     Ke - (2x2) matrice de raideur diffusion
%     Me - (2x2) matrice de masse
%
%   Voir aussi : elastic_element_stiff, solve_transient_cyclic

L = x2 - x1;
if L <= 0
    error('diffusion_element_stiff:negativeLength', ...
          'Longueur d''élément négative ou nulle : L = %e', L);
end

Ke = (D / L) * [1, -1; -1, 1];
Me = (L / 6) * [2, 1; 1, 2];
end
