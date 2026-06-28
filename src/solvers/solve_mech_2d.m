function sol = solve_mech_2d(P, model, BC, opts)
%SOLVE_MECH_2D  2D plane-strain chemo-mechanical solver (FULL or EFFECTIVE).
%
%   sol = solve_mech_2d(P, model, BC, opts) solves the static plane-strain
%   elasticity problem (optionally with a uniform compositional eigenstrain in a
%   chosen region) for the trilayer. model = 'full' (SEI meshed) or 'eff'
%   (SEI replaced by an anisotropic interface with normal stiffness K_n and
%   tangential stiffness K_t).
%
%   BC fields (each a struct array with .nodes and .val):
%     BC.ux : prescribed x-displacement on node list   (optional)
%     BC.uy : prescribed y-displacement on node list   (optional)
%   opts : mesh options (see build_*_mesh_2d), plus
%     opts.Kn, opts.Kt : effective interface stiffnesses [Pa/m] (EFF model)
%     opts.swell_region : region id with eigenstrain (optional)
%     opts.dc : concentration change for the eigenstrain (optional)
%
%   Output sol: .nodes, .u (Nn x 2), .m (mesh), .K (assembled, pre-BC), .F.
%
%   See also: elastic_q4_2d, build_trilayer_mesh_2d, build_eff_mesh_2d

if strcmpi(model,'full')
    m = build_trilayer_mesh_2d(P, opts);
else
    m = build_eff_mesh_2d(P, opts);
end
Nn = size(m.nodes,1); ndof = 2*Nn;

swell_reg = getf(opts,'swell_region', 0);
dc_val    = getf(opts,'dc', 0);

% ----- assemble bulk elasticity -----
I = []; J = []; V = []; F = zeros(ndof,1);
for e = 1:size(m.elems,1)
    en = m.elems(e,:); coords = m.nodes(en,:);
    [E,nu,beta] = layer_props(P, m.eregion(e));
    dce = (m.eregion(e)==swell_reg) * dc_val;
    [Ke,fe] = elastic_q4_2d(coords, E, nu, beta, dce);
    dofs = reshape([2*en-1; 2*en], 1, []);
    [aa,bb] = meshgrid(dofs,dofs);
    I = [I; bb(:)]; J = [J; aa(:)]; V = [V; Ke(:)]; %#ok<AGROW>
    F(dofs) = F(dofs) + fe;
end
K = sparse(I, J, V, ndof, ndof);

% ----- effective interface (anisotropic springs) -----
if strcmpi(model,'eff')
    Kn = opts.Kn; Kt = opts.Kt;
    beff = getf(opts,'beta_eff',0);
    xb = m.iface_x;
    for k = 1:numel(m.iface_bot)
        nb = m.iface_bot(k); nt = m.iface_top(k);
        % tributary length along the interface
        if k==1,                Ltr = 0.5*(xb(2)-xb(1));
        elseif k==numel(xb),    Ltr = 0.5*(xb(end)-xb(end-1));
        else,                   Ltr = 0.5*(xb(k+1)-xb(k-1));
        end
        % tangential spring (x dofs), normal spring (y dofs)
        for d = [1 2]
            if d==1, kk = Kt*Ltr; db = 2*nb-1; dt = 2*nt-1;   % x
            else,    kk = Kn*Ltr; db = 2*nb;   dt = 2*nt;     % y
            end
            K(db,db)=K(db,db)+kk; K(dt,dt)=K(dt,dt)+kk;
            K(db,dt)=K(db,dt)-kk; K(dt,db)=K(dt,db)-kk;
        end
        % chemo-mechanical swelling force on the normal dofs
        if beff~=0 && dc_val~=0
            f = Kn*Ltr*beff*dc_val;
            F(2*nb)=F(2*nb)-f; F(2*nt)=F(2*nt)+f;
        end
    end
end

% ----- Dirichlet BCs (val may be a scalar or a per-node vector) -----
fixed = []; val = [];
if isfield(BC,'ux') && ~isempty(BC.ux)
    nn = BC.ux.nodes(:)'; vv = BC.ux.val; if isscalar(vv), vv = vv*ones(1,numel(nn)); end
    fixed=[fixed, 2*nn-1]; val=[val, vv(:)'];
end
if isfield(BC,'uy') && ~isempty(BC.uy)
    nn = BC.uy.nodes(:)'; vv = BC.uy.val; if isscalar(vv), vv = vv*ones(1,numel(nn)); end
    fixed=[fixed, 2*nn];   val=[val, vv(:)'];
end
fixed = fixed(:)'; val = val(:)';
[fixed, ia] = unique(fixed,'last'); val = val(ia);

% ----- optional periodic-x reduction (master-slave: u(right)=u(left)) -----
if getf(opts,'periodic_x',false)
    L = m.left(:)'; R = m.right(:)';
    if numel(L) ~= numel(R)
        error('solve_mech_2d:periodic','left/right node counts differ.');
    end
    sl = reshape([2*R-1; 2*R], 1, []);     % slave dofs (right)
    ms = reshape([2*L-1; 2*L], 1, []);     % master dofs (left)
    rep = 1:ndof; rep(sl) = ms;            % representative dof
    keep = setdiff(1:ndof, sl);
    redidx = zeros(1,ndof); redidx(keep) = 1:numel(keep);
    T = sparse(1:ndof, redidx(rep), 1, ndof, numel(keep));
    Kr = T'*K*T; Fr = T'*F;
    fr = redidx(rep(fixed));
    [fr, ib] = unique(fr,'last'); vr = val(ib);
    ndr = numel(keep);
else
    T = speye(ndof); Kr = K; Fr = F; fr = fixed; vr = val; ndr = ndof;
end

free = setdiff((1:ndr)', fr(:));
ur = zeros(ndr,1); ur(fr) = vr;
Fr = Fr - Kr*ur;
ur(free) = Kr(free,free) \ Fr(free);
u = T*ur;

sol.nodes = m.nodes;
sol.u = [u(1:2:end), u(2:2:end)];
sol.m = m; sol.K = K; sol.uvec = u;
end

% --------------------------------------------------------------------------
function [E,nu,beta] = layer_props(P, reg)
switch reg
    case 1, E=P.E1;    nu=P.nu1;   beta=P.beta1;     % electrode
    case 2, E=P.E_sei; nu=P.nu_sei;beta=P.beta_sei;  % SEI (inner sublayer / homogeneous)
    case 3, E=P.E2;    nu=P.nu2;   beta=P.beta2;     % electrolyte
    case 4                                            % SEI outer sublayer (bilayer)
        E=getf2(P,'E_sei2',P.E_sei); nu=getf2(P,'nu_sei2',P.nu_sei); beta=getf2(P,'beta_sei2',P.beta_sei);
end
end

function v = getf2(s,f,d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function v = getf(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
