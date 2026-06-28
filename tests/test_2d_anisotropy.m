function out = test_2d_anisotropy()
%TEST_2D_ANISOTROPY  2D plane-strain validation of the normal/tangential
%   interface anisotropy K_n/K_t -- the key result the 1D model cannot test.
%
%   Symmetric trilayer  bulk | SEI | bulk  (both bulks = electrode properties,
%   so the SEI is the only special layer), plane strain.
%
%   (A) NORMAL (roller sides, affine uy = eps0*y) -> K_n = sigma_yy/[uy] vs (lam+2mu)/eps.
%   (B) SHEAR  (affine KUBC u=(gamma*y,0))        -> K_t = sigma_xy/[ux] vs mu/eps.
%   Anisotropy K_n/K_t is compared to 2(1-nu)/(1-2nu); an ISOTROPIC interface
%   (K_t=K_n) is shown to fail in shear, proving the anisotropy is necessary.

addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src')));

fprintf('\n===============================================================\n');
fprintf('  TEST 2D : interface anisotropy  K_n / K_t  (plane strain)\n');
fprintf('===============================================================\n');

P = config_parameters('intermediate');
P.E2 = P.E1; P.nu2 = P.nu1; P.beta2 = 0;          % symmetric bulk
lam = P.E_sei*P.nu_sei/((1+P.nu_sei)*(1-2*P.nu_sei));
mu  = P.E_sei/(2*(1+P.nu_sei));
Kn_th = (lam+2*mu)/P.eps;  Kt_th = mu/P.eps;  ratio_th = (lam+2*mu)/mu;

opts = struct('Lx',12*P.eps,'Le',max(20*P.eps,6e-7),'nx',24,'nye',10,'nys',6,'nyel',10);
opts.Kn = Kn_th; opts.Kt = Kt_th; opts.beta_eff = 0;
Le = opts.Le; ep = P.eps; Lx = opts.Lx;

mf = build_trilayer_mesh_2d(P,opts);  me = build_eff_mesh_2d(P,opts);

% ---------------- (A) NORMAL : rollers + affine uy = eps0*y ----------------
eps0 = 1e-3;
BCf = struct();
BCf.ux = struct('nodes',[mf.left,mf.right],'val',0);
BCf.uy = struct('nodes',[mf.bottom,mf.top], ...
   'val',[zeros(1,numel(mf.bottom)), eps0*mf.nodes(mf.top,2)']);
solF = solve_mech_2d(P,'full',BCf,opts);

BCe = struct();
BCe.ux = struct('nodes',[me.left,me.right],'val',0);
BCe.uy = struct('nodes',[me.bottom,me.top], ...
   'val',[zeros(1,numel(me.bottom)), eps0*me.nodes(me.top,2)']);
solE = solve_mech_2d(P,'eff',BCe,opts);

% K_n from FULL : sigma_yy = reaction/Lx ; normal jump across SEI at mid-x
Rn = solF.K*solF.uvec;
sig_yy = sum(Rn(2*mf.top))/Lx;
Dn_full = sample(solF,mf,Lx/2,Le+ep,2) - sample(solF,mf,Lx/2,Le,2);
Kn_meas = sig_yy/Dn_full;
% normal jump reproduced by the EFF interface (doubled nodes at y=Le)
Dn_eff  = ifjump(solE,Lx/2,2);
fprintf('\n[A] NORMAL   K_n meas = %.3e   theory (lam+2mu)/eps = %.3e   (ratio %.4f)\n',...
        Kn_meas,Kn_th,Kn_meas/Kn_th);
fprintf('             normal jump [uy]:  FULL = %.3e   EFF(anis) = %.3e   (ratio %.3f)\n',...
        Dn_full,Dn_eff,Dn_eff/Dn_full);

% ---------------- (B) SHEAR : periodic-x simple shear ----------------
gamma = 1e-3;
optsP = opts; optsP.periodic_x = true;              % uniform-in-x -> clean shear
optsPi = optsP; optsPi.Kt = Kn_th;                  % isotropic control (K_t=K_n)
solFs   = solve_mech_2d(P,'full',shearBC(mf,gamma),optsP);
solEs   = solve_mech_2d(P,'eff', shearBC(me,gamma),optsP);
solEiso = solve_mech_2d(P,'eff', shearBC(me,gamma),optsPi);

% sigma_xy = reaction_x(top)/Lx  (clean under periodicity); jump across SEI
Rs = solFs.K*solFs.uvec;
sig_xy = sum(Rs(2*mf.top-1))/Lx;
Dt_full = sample(solFs,mf,Lx/2,Le+ep,1) - sample(solFs,mf,Lx/2,Le,1);
Kt_meas = sig_xy/Dt_full;
Dt_anis = ifjump(solEs,  Lx/2,1);
Dt_iso  = ifjump(solEiso,Lx/2,1);
fprintf('\n[B] SHEAR    K_t meas = %.3e   theory mu/eps        = %.3e   (ratio %.4f)\n',...
        Kt_meas,Kt_th,Kt_meas/Kt_th);
fprintf('             tangential jump [ux]:  FULL = %.3e\n', Dt_full);
fprintf('               EFF anisotropic (K_t=mu/eps)   = %.3e   (ratio %.3f, should ~1)\n',...
        Dt_anis, Dt_anis/Dt_full);
fprintf('               EFF ISOTROPIC   (K_t=K_n)      = %.3e   (ratio %.3f, should ~%.2f)\n',...
        Dt_iso, Dt_iso/Dt_full, Kt_th/Kn_th);

fprintf('\n[C] ANISOTROPY  K_n/K_t :  measured = %.3f   theory 2(1-nu)/(1-2nu) = %.3f\n',...
        Kn_meas/Kt_meas, ratio_th);
fprintf('DONE.\n');

out = struct('Kn_meas',Kn_meas,'Kn_th',Kn_th,'Kt_meas',Kt_meas,'Kt_th',Kt_th,...
   'ratio_meas',Kn_meas/Kt_meas,'ratio_th',ratio_th,...
   'Dn_full',Dn_full,'Dn_eff',Dn_eff,'Dt_full',Dt_full,'Dt_anis',Dt_anis,'Dt_iso',Dt_iso);
end

function dj = ifjump(sol, xq, comp)
% interface jump in the EFF mesh: u(top side) - u(bottom side) at x ~ xq
m=sol.m; xb=m.iface_x; [~,k]=min(abs(xb-xq));
dj = sol.u(m.iface_top(k),comp) - sol.u(m.iface_bot(k),comp);
end

% ---------- helpers ----------
function BC = shearBC(m, gamma)
% simple shear with periodic x: bottom fixed, top ux = gamma*Htop, uy=0
Htop = max(m.nodes(:,2));
BC.ux = struct('nodes',[m.bottom,m.top], ...
   'val',[zeros(1,numel(m.bottom)), gamma*Htop*ones(1,numel(m.top))]);
BC.uy = struct('nodes',[m.bottom,m.top],'val',0);
end
function v = sample(sol,m,xq,yq,comp)
d=(m.nodes(:,1)-xq).^2+(m.nodes(:,2)-yq).^2; [~,k]=min(d); v=sol.u(k,comp);
end
function e = electrode_L2(solA,solB,comp,Le)
mA=solA.m; idx=find(mA.nodes(:,2)<=Le+1e-12);
ua=solA.u(idx,comp); ub=zeros(size(ua));
for q=1:numel(idx)
    d=(solB.m.nodes(:,1)-mA.nodes(idx(q),1)).^2+(solB.m.nodes(:,2)-mA.nodes(idx(q),2)).^2;
    [~,k]=min(d); ub(q)=solB.u(k,comp);
end
num=sqrt(sum((ua-ub).^2)); den=sqrt(sum(ua.^2)); if den<1e-30,den=1;end; e=num/den;
end
