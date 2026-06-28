function m = build_trilayer_mesh_2d(P, opts)
%BUILD_TRILAYER_MESH_2D  Structured Q4 mesh of the FULL 2D trilayer.
%
%   m = build_trilayer_mesh_2d(P, opts) builds a structured quadrilateral mesh of
%   the electrode / SEI / electrolyte trilayer in plane strain. The interface
%   normal is the +y direction; x is the tangential direction.
%
%     electrolyte  y in [Le+eps, Le+eps+Le]   (region 3)
%     SEI          y in [Le, Le+eps]           (region 2)   <- explicitly meshed
%     electrode    y in [0, Le]                (region 1)
%
%   opts fields (with defaults): Lx (tangential length), Le (bulk height),
%   nx, nye, nys, nyel (element counts).
%
%   Output struct m:
%     .nodes (Nn x 2), .elems (Ne x 4, CCW), .eregion (Ne x 1)
%     .bottom, .top, .left, .right : boundary node index lists
%     .Lx, .Le, .eps
%
%   See also: build_eff_mesh_2d, elastic_q4_2d

eps = P.eps;
Lx  = getf(opts,'Lx', 10*eps);
Le  = getf(opts,'Le', max(20*eps, 5e-7));
nx  = getf(opts,'nx', 24);
nye = getf(opts,'nye', 10);
nys = getf(opts,'nys', 6);
nyel= getf(opts,'nyel',10);
h1f = getf(opts,'h1_frac', 0);          % >0 : bilayer SEI, inner sublayer = h1f*eps
h1  = h1f*eps;

xv = linspace(0, Lx, nx+1);
if h1f > 0   % bilayer: mesh the inner/outer sublayer boundary at Le+h1
    n1 = max(2,round(nys*h1f)); n2 = max(2,nys-n1);
    seiy = [linspace(Le,Le+h1,n1+1), linspace(Le+h1,Le+eps,n2+1)];
else
    seiy = linspace(Le,Le+eps,nys+1);
end
yv = unique([linspace(0,Le,nye+1), seiy, linspace(Le+eps,Le+eps+Le,nyel+1)]);
nxn = numel(xv); nyn = numel(yv);

% Nodes (tensor product), id(i,j) = (j-1)*nxn + i
[X,Y] = meshgrid(xv, yv);          % nyn x nxn
nodes = [reshape(X',[],1), reshape(Y',[],1)];   % ordered i fastest

id = @(i,j) (j-1)*nxn + i;

% Elements + region
elems = zeros((nxn-1)*(nyn-1),4);
ereg  = zeros((nxn-1)*(nyn-1),1);
e = 0;
for j = 1:nyn-1
    yc = 0.5*(yv(j)+yv(j+1));
    if yc < Le - 1e-15
        reg = 1;                                   % electrode
    elseif yc < Le+eps - 1e-15
        if h1f > 0 && yc >= Le+h1 - 1e-15
            reg = 4;                               % SEI outer sublayer (organic)
        else
            reg = 2;                               % SEI inner sublayer (inorganic) / homogeneous SEI
        end
    else
        reg = 3;                                   % electrolyte
    end
    for i = 1:nxn-1
        e = e + 1;
        elems(e,:) = [id(i,j), id(i+1,j), id(i+1,j+1), id(i,j+1)];
        ereg(e)    = reg;
    end
end

m.nodes = nodes; m.elems = elems; m.eregion = ereg;
m.bottom = id((1:nxn),1);
m.top    = id((1:nxn),nyn);
m.left   = id(1,(1:nyn));
m.right  = id(nxn,(1:nyn));
m.Lx = Lx; m.Le = Le; m.eps = eps;
m.nxn = nxn; m.nyn = nyn; m.xv = xv; m.yv = yv;
end

function v = getf(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
