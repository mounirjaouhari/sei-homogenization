function m = build_eff_mesh_2d(P, opts)
%BUILD_EFF_MESH_2D  Structured Q4 mesh of the EFFECTIVE 2D model (SEI -> interface).
%
%   m = build_eff_mesh_2d(P, opts) builds electrode and electrolyte bulk meshes
%   joined by a zero-thickness interface line at y = Le. The interface nodes are
%   DOUBLED: for each tangential position x_i there is a bottom node (electrode
%   side) and a top node (electrolyte side) at the same coordinate (x_i, Le).
%
%   Output struct m (same bulk layout as build_trilayer_mesh_2d, SEI removed):
%     .nodes, .elems, .eregion, .bottom, .top, .left, .right
%     .iface_bot (1 x nxn), .iface_top (1 x nxn) : paired interface node indices
%     .iface_x   (1 x nxn) : tangential coordinates of the interface nodes
%
%   See also: build_trilayer_mesh_2d, interface_element_2d

eps = P.eps;
Lx  = getf(opts,'Lx', 10*eps);
Le  = getf(opts,'Le', max(20*eps, 5e-7));
nx  = getf(opts,'nx', 24);
nye = getf(opts,'nye', 10);
nyel= getf(opts,'nyel',10);

xv = linspace(0, Lx, nx+1);
nxn = numel(xv);

% Electrode nodes (y in [0,Le]) then DOUBLED interface row, then electrolyte.
ye  = linspace(0, Le, nye+1);          % includes y=Le (last)
yel = linspace(Le, Le+Le, nyel+1);     % includes y=Le (first)

% Electrode block: rows 1..nye+1 (y = ye), nxn columns
% Electrolyte block: rows starting after, y = yel
% The shared y=Le is represented TWICE (electrode top row, electrolyte bottom row).
Xe = repmat(xv', numel(ye),1);
Ye = reshape(repmat(ye, nxn,1),[],1);
Xl = repmat(xv', numel(yel),1);
Yl = reshape(repmat(yel,nxn,1),[],1);
nodes = [ [Xe,Ye]; [Xl,Yl] ];

nE = numel(ye)*nxn;                     % number of electrode nodes
ide = @(i,j) (j-1)*nxn + i;             % electrode node id (j=1..nye+1)
idl = @(i,j) nE + (j-1)*nxn + i;        % electrolyte node id (j=1..nyel+1)

% Elements
elems = []; ereg = [];
for j = 1:numel(ye)-1
    for i = 1:nxn-1
        elems(end+1,:) = [ide(i,j), ide(i+1,j), ide(i+1,j+1), ide(i,j+1)]; %#ok<AGROW>
        ereg(end+1,1)  = 1; %#ok<AGROW>
    end
end
for j = 1:numel(yel)-1
    for i = 1:nxn-1
        elems(end+1,:) = [idl(i,j), idl(i+1,j), idl(i+1,j+1), idl(i,j+1)]; %#ok<AGROW>
        ereg(end+1,1)  = 3; %#ok<AGROW>
    end
end

% Interface node pairs (electrode top row = ye end; electrolyte bottom row = yel 1)
iface_bot = ide((1:nxn), numel(ye));    % electrode side
iface_top = idl((1:nxn), 1);            % electrolyte side

m.nodes = nodes; m.elems = elems; m.eregion = ereg;
m.iface_bot = iface_bot; m.iface_top = iface_top; m.iface_x = xv;
m.bottom = ide((1:nxn),1);
m.top    = idl((1:nxn),numel(yel));
m.left   = [ide(1,(1:numel(ye))), idl(1,(1:numel(yel)))];
m.right  = [ide(nxn,(1:numel(ye))), idl(nxn,(1:numel(yel)))];
m.Lx = Lx; m.Le = Le; m.eps = eps; m.nxn = nxn; m.xv = xv;
end

function v = getf(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
