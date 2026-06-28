function [nodes, elems, regions] = build_eff_mesh(P)
%BUILD_EFF_MESH  Maillage 1D du modèle EFFECTIF (SEI remplacée par interface)
%
%   [nodes, elems, regions] = build_eff_mesh(P) génère un maillage 1D
%   pour le modèle effectif :
%
%      électrode | interface (zéro épaisseur) | électrolyte
%      <---L--->                              <---L--->
%
%   La SEI est remplacée par une interface d'épaisseur nulle située au
%   nœud d'indice idx_interface, portant les conditions de saut dérivées
%   dans Theorem 3.1.
%
%   Sorties :
%     nodes   - (Nn x 1) coordonnées des nœuds (m)
%     elems   - (Ne x 2) connectivité élémentaire
%     regions - structure avec champs :
%       .elems                - (Ne x 1) région par élément
%       .idx_interface_left   - dernier nœud côté électrode
%       .idx_interface_right  - premier nœud côté électrolyte
%       .idx_left, .idx_right - bords
%       .n_nodes, .n_elems
%
%   Voir aussi : build_trilayer_mesh, assemble_eff_model

% =====================================================================
%  Maillage de l'électrode (région 1)
% =====================================================================
L1 = P.L;
x1 = linspace(0, L1, P.n_elem_electrode + 1)';

% Maillage de l'électrolyte (région 3) — la SEI est absente
L2 = P.L;
x2 = linspace(L1, L1 + L2, P.n_elem_electrolyte + 1)';

% Concaténation : le nœud L1 est partagé entre les deux maillages
% mais correspond à DEUX nœuds distincts dans le modèle effectif
% (nœud gauche de l'interface et nœud droit de l'interface)
nodes = [x1; x2(2:end)];
Nn = length(nodes);

% Indice du nœud interface côté électrode (dernier nœud de l'électrode)
idx_left_iface  = P.n_elem_electrode + 1;
% Indice du nœud interface côté électrolyte (premier nœud de l'électrolyte)
idx_right_iface = idx_left_iface;  % même coordonnée mais sera dédoublée

% En fait, pour le modèle effectif, on crée DEUX nœuds au point z=L1 :
%   - nœud idx_left_iface  : extrémité droite de l'électrode
%   - nœud idx_right_iface : extrémité gauche de l'électrolyte
% On reconstruit nodes en intercalant le nœud dédoublé :
nodes = [x1(1:end-1);          % électrode sans son dernier nœud
         x1(end);              % nœud gauche interface (électrode)
         x1(end);              % nœud droit interface (électrolyte)
         x2(2:end)];           % électrolyte sans son premier nœud
Nn = length(nodes);

% Le tableau nodes = [x1(1:end-1); x1(end); x1(end); x2(2:end)] :
%   indices 1..n_elem_electrode      -> électrode (0 .. L1-dx)
%   index   n_elem_electrode+1       -> nœud L1 GAUCHE (fin électrode)
%   index   n_elem_electrode+2       -> nœud L1 DROIT  (début électrolyte)
%   indices n_elem_electrode+3..Nn   -> électrolyte
idx_left_iface  = P.n_elem_electrode + 1;    % nœud L1 côté électrode
idx_right_iface = P.n_elem_electrode + 2;    % nœud L1 côté électrolyte (dédoublé)
% Note : nodes(idx_left_iface) == nodes(idx_right_iface) == L1

% =====================================================================
%  Connectivité élémentaire (éléments linéaires 2-nœuds)
% =====================================================================
elems = [];
regions_list = [];

% Électrode : éléments 1 à n_elem_electrode
for i = 1:P.n_elem_electrode
    elems = [elems; i, i+1]; %#ok<AGROW>
    regions_list = [regions_list; 1]; %#ok<AGROW>
end

% Électrolyte : éléments n_elem_electrode+1 à n_elem_electrode+n_elem_electrolyte
for i = 1:P.n_elem_electrolyte
    elems = [elems; idx_right_iface + i - 1, idx_right_iface + i]; %#ok<AGROW>
    regions_list = [regions_list; 3]; %#ok<AGROW>
end

% =====================================================================
%  Structure des régions et indices
% =====================================================================
regions = struct(...
    'elems', regions_list, ...
    'idx_interface_left', idx_left_iface, ...
    'idx_interface_right', idx_right_iface, ...
    'idx_left', 1, ...
    'idx_right', Nn, ...
    'n_nodes', Nn, ...
    'n_elems', size(elems, 1));

assert(size(elems, 1) == P.n_elem_electrode + P.n_elem_electrolyte, ...
       'build_eff_mesh:inconsistent', 'Incohérence dans le nombre d''éléments');
end
