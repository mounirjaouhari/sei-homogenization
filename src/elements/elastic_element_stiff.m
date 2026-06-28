function [Ke, fe] = elastic_element_stiff(x1, x2, C1d, beta, c_elem)
%ELASTIC_ELEMENT_STIFF  Matrice de raideur élémentaire pour élasticité 1D
%                       couplée à la diffusion (eigenstrain chimique).
%
%   [Ke, fe] = elastic_element_stiff(x1, x2, C1d, beta, c_elem) calcule
%   la matrice de raideur (2x2) et le vecteur forces (2x1) pour un élément
%   fini linéaire 1D à deux nœuds, avec :
%
%     - Module élastique 1D : C1d = lambda + 2*mu (Pa)
%     - Coefficient de Vegard : beta (m^3/mol)
%     - Concentration aux nœuds : c_elem = [c1; c2] (mol/m^3)
%
%   La formulation faible est :
%
%     int_Omega sigma : eps(v) dV = int_Omega f v dV
%
%   avec sigma = C1d * (du/dx - beta * (c - c0))
%
%   Discrétisation : u = N1*u1 + N2*u2, c = N1*c1 + N2*c2
%   avec N1 = (x2-x)/L, N2 = (x-x1)/L, L = x2-x1
%
%   La matrice élémentaire se décompose :
%     Ke = K_elastic + K_chemo
%   où K_elastic = (C1d/L) * [[1, -1]; [-1, 1]]
%   et le couplage chémo-mécanique est dans fe (force équivalente).
%
%   Sorties :
%     Ke - (2x2) matrice de raideur
%     fe - (2x1) vecteur force équivalent (contribution de l'eigenstrain)
%
%   Voir aussi : diffusion_element_stiff, assemble_full_model

L = x2 - x1;
if L <= 0
    error('elastic_element_stiff:negativeLength', ...
          'Longueur d''élément négative ou nulle : L = %e', L);
end

% Matrice de raideur élastique pure (1D, linéaire)
Ke = (C1d / L) * [1, -1; -1, 1];

% Force équivalente due à l'eigenstrain chimique :
%   sigma = C1d * (du/dx - beta * (c - c0))
%   => terme de force : int_Omega C1d * beta * (c - c0) * dN/dx dV
%   Pour c = N1*c1 + N2*c2 :
%     int_beta = C1d * beta * [ (c1-c0)*int(dN/dx*N1) + (c2-c0)*int(dN/dx*N2) ]
%   int_0^L dN1/dx * N1 dx = -1/2,  int_0^L dN1/dx * N2 dx = 1/2
%   int_0^L dN2/dx * N1 dx = -1/2, int_0^L dN2/dx * N2 dx = 1/2
% Donc fe = C1d * beta * [-(c1+c2)/2 - c0*(-1+1)/2 ;  (c1+c2)/2 - c0*0]
%         = C1d * beta * (c_avg - c0) * [-1 ; 1] / 1
% avec c_avg = (c1 + c2)/2

c_avg = 0.5 * (c_elem(1) + c_elem(2));
fe = C1d * beta * (c_avg - 0) * [-1; 1];  % c0 = 0 par défaut

% Note : c0 est supposé nul ici (référence à l'état vierge)
end
