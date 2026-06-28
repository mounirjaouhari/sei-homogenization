function out = test_2d_convergence()
%TEST_2D_CONVERGENCE  2D convergence of the effective interface model to the
%   resolved trilayer, in BOTH normal and shear (tangential) modes, as eps->0.
%
%   For a range of eps/L, the resolved (FULL) and effective (EFF, anisotropic
%   K_n,K_t) models are solved under (i) normal uniaxial strain and (ii) affine
%   shear KUBC; the relative error on the interface jump is fitted vs eps.

addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src')));
fprintf('\n===============================================================\n');
fprintf('  TEST 2D : convergence (normal + shear) of the interface model\n');
fprintf('===============================================================\n');

P = config_parameters('intermediate');
P.E2=P.E1; P.nu2=P.nu1; P.beta2=0;
lam=P.E_sei*P.nu_sei/((1+P.nu_sei)*(1-2*P.nu_sei)); mu=P.E_sei/(2*(1+P.nu_sei));

Le = 6e-7;  Lfac = 12;
eps_list = Le*[0.02 0.04 0.08 0.16];      % eps/Le
errN = zeros(size(eps_list)); errT = errN;

for i=1:numel(eps_list)
    P.eps = eps_list(i); ep=P.eps;
    opts = struct('Lx',Lfac*ep,'Le',Le,'nx',20,'nye',10,'nys',6,'nyel',10);
    opts.Kn=(lam+2*mu)/ep; opts.Kt=mu/ep; opts.beta_eff=0;
    mf=build_trilayer_mesh_2d(P,opts); me=build_eff_mesh_2d(P,opts); Lx=opts.Lx;
    % normal
    eps0=1e-3;
    bcF=struct('ux',struct('nodes',[mf.left,mf.right],'val',0), ...
        'uy',struct('nodes',[mf.bottom,mf.top],'val',[zeros(1,numel(mf.bottom)),eps0*mf.nodes(mf.top,2)']));
    bcE=struct('ux',struct('nodes',[me.left,me.right],'val',0), ...
        'uy',struct('nodes',[me.bottom,me.top],'val',[zeros(1,numel(me.bottom)),eps0*me.nodes(me.top,2)']));
    sF=solve_mech_2d(P,'full',bcF,opts); sE=solve_mech_2d(P,'eff',bcE,opts);
    DnF=samp(sF,Lx/2,Le+ep,2)-samp(sF,Lx/2,Le,2); DnE=ifj(sE,Lx/2,2);
    errN(i)=abs(DnE-DnF)/abs(DnF);
    % shear (periodic x -> clean simple shear)
    g=1e-3; optsP=opts; optsP.periodic_x=true;
    sFs=solve_mech_2d(P,'full',shr(mf,g),optsP); sEs=solve_mech_2d(P,'eff',shr(me,g),optsP);
    DtF=samp(sFs,Lx/2,Le+ep,1)-samp(sFs,Lx/2,Le,1); DtE=ifj(sEs,Lx/2,1);
    errT(i)=abs(DtE-DtF)/abs(DtF);
    fprintf('  eps/Le=%.3f : normal jump err=%.3e , shear jump err=%.3e\n', ...
            ep/Le, errN(i), errT(i));
end
pN=polyfit(log(eps_list/Le),log(errN),1); pT=polyfit(log(eps_list/Le),log(errT),1);
fprintf('\n  Convergence rate (jump error vs eps):  normal = %.2f , shear = %.2f\n',pN(1),pT(1));
out=struct('eps_over_Le',eps_list/Le,'errN',errN,'errT',errT,'rateN',pN(1),'rateT',pT(1));
fprintf('DONE.\n');
end
function v=samp(sol,xq,yq,c), m=sol.m; d=(m.nodes(:,1)-xq).^2+(m.nodes(:,2)-yq).^2;[~,k]=min(d);v=sol.u(k,c); end
function dj=ifj(sol,xq,c), m=sol.m;[~,k]=min(abs(m.iface_x-xq)); dj=sol.u(m.iface_top(k),c)-sol.u(m.iface_bot(k),c); end
function BC=shr(m,g), Htop=max(m.nodes(:,2)); BC.ux=struct('nodes',[m.bottom,m.top],'val',[zeros(1,numel(m.bottom)),g*Htop*ones(1,numel(m.top))]); BC.uy=struct('nodes',[m.bottom,m.top],'val',0); end
