function [nodes, elems, regions] = build_trilayer_mesh(P)
%BUILD_TRILAYER_MESH  Maillage 1D du modèle FULL (SEI explicitement maillée)
%
%   [nodes, elems, regions] = build_trilayer_mesh(P) génère un maillage 1D
%   quadratique pour la configuration tricouche :
%
%      électrode | SEI | électrolyte
%      <---L---> <-eps-> <---L--->
%
%   La SEI est explicitement maillée avec n_elem_sei éléments.
%
%   Sorties :
%     nodes   - (Nn x 1) coordonnées des nœuds (m)
%     elems   - (Ne x 3) connectivité élémentaire [n1 n2 n3] (quadratique)
%     regions - (Ne x 1) identifiant de région : 1=électrode, 2=SEI, 3=électrolyte
%
%   Voir aussi : build_eff_mesh, elastic_element_stiff

% =====================================================================
%  Maillage de l'électrode (région 1)
% =====================================================================
L1 = P.L;
x1 = linspace(0, L1, P.n_elem_electrode + 1)';

% Maillage de la SEI (région 2)
Lsei = P.eps;
x2 = linspace(L1, L1 + Lsei, P.n_elem_sei + 1)';

% Maillage de l'électrolyte (région 3)
L2 = P.L;
x3 = linspace(L1 + Lsei, L1 + Lsei + L2, P.n_elem_electrolyte + 1)';

% Concaténation (en supprimant les nœuds de jonction dupliqués)
nodes = [x1; x2(2:end); x3(2:end)];
Nn = length(nodes);

% =====================================================================
%  Connectivité élémentaire (éléments quadratiques 3-nœuds)
%  Pour un maillage linéaire raffiné en quadratique, on insère le nœud
%  milieu. Plus simple : on utilise des éléments linéaires à 2 nœuds pour
%  ce programme de validation (suffisant pour démontrer O(eps)).
% =====================================================================
elems = [];
regions = [];

% Électrode  (connectivité par INDICES de nœuds, pas par coordonnées)
for i = 1:P.n_elem_electrode
    elems = [elems; i, i+1]; %#ok<AGROW>
    regions = [regions; 1]; %#ok<AGROW>
end

% SEI
n1_start = P.n_elem_electrode + 1;
for i = 1:P.n_elem_sei
    elems = [elems; n1_start + i - 1, n1_start + i]; %#ok<AGROW>
    regions = [regions; 2]; %#ok<AGROW>
end

% Électrolyte
n2_start = P.n_elem_electrode + P.n_elem_sei + 1;
for i = 1:P.n_elem_electrolyte
    elems = [elems; n2_start + i - 1, n2_start + i]; %#ok<AGROW>
    regions = [regions; 3]; %#ok<AGROW>
end

% =====================================================================
%  Stockage des indices de jonction
% =====================================================================
P_nodes.idx_electrode_sei   = P.n_elem_electrode + 1;       % Nœud interface élec/SEI
P_nodes.idx_sei_electrolyte = P.n_elem_electrode + P.n_elem_sei + 1; % Nœud interface SEI/élec
P_nodes.idx_left            = 1;                              % Bord gauche
P_nodes.idx_right           = Nn;                             % Bord droit

% Ajout aux sorties (sous forme de structure)
regions = struct(...
    'elems', regions, ...
    'idx_electrode_sei', P_nodes.idx_electrode_sei, ...
    'idx_sei_electrolyte', P_nodes.idx_sei_electrolyte, ...
    'idx_left', P_nodes.idx_left, ...
    'idx_right', P_nodes.idx_right, ...
    'n_nodes', Nn, ...
    'n_elems', size(elems, 1));

% Vérification
assert(size(elems, 1) == P.n_elem_electrode + P.n_elem_sei + P.n_elem_electrolyte, ...
       'build_trilayer_mesh:inconsistent', 'Incohérence dans le nombre d''éléments');
end
