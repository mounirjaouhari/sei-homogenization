function out = test_2d_bilayer()
%TEST_2D_BILAYER  2D plane-strain validation of the BILAYER interface law.
%
%   The real SEI is a bilayer: a stiff, Li-blocking inorganic inner sublayer
%   (LiF/Li2CO3) and a soft, Li-conducting organic outer sublayer. Treating the
%   two sublayers as springs in series through the thickness gives the effective
%   interface compliances
%       A_n = sum_i h_i/(lam_i+2mu_i),   A_t = sum_i h_i/mu_i,   K = 1/A,
%       R_eff = sum_i h_i/D_i .
%   The NON-ADDITIVE consequence: the mechanical stiffness is controlled by the
%   SOFT (organic) layer, while the diffusion resistance is controlled by the
%   BLOCKING (inorganic) layer -- the two effective properties are governed by
%   DIFFERENT sublayers. This test validates K_n, K_t against the series formula
%   in plane strain (resolved bilayer FULL vs effective interface).

addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src')));
fprintf('\n===============================================================\n');
fprintf('  TEST 2D : BILAYER interface (inorganic inner + organic outer)\n');
fprintf('===============================================================\n');

P = config_parameters('intermediate');
P.E2=P.E1; P.nu2=P.nu1; P.beta2=0;                       % symmetric bulk
% --- bilayer SEI sublayers (inner inorganic, outer organic), each half-thickness ---
P.E_sei  = 50e9;  P.nu_sei  = 0.24; P.beta_sei  = 2e-3;  % inner : LiF-like (stiff, blocking)
P.E_sei2 = 2e9;   P.nu_sei2 = 0.35; P.beta_sei2 = 8e-3;  % outer : organic (soft, conducting)
P.D_sei  = 1e-17; P.D_sei2 = 1e-14;                      % inner blocks Li, outer conducts
h1f = 0.5;  ep = P.eps;  h1 = h1f*ep;  h2 = ep-h1;

lame = @(E,nu) deal(E*nu/((1+nu)*(1-2*nu)), E/(2*(1+nu)));
[l1,m1]=lame(P.E_sei ,P.nu_sei );  [l2,m2]=lame(P.E_sei2,P.nu_sei2);
% series (harmonic) effective stiffnesses
Kn_th = 1/( h1/(l1+2*m1) + h2/(l2+2*m2) );
Kt_th = 1/( h1/m1        + h2/m2        );
% reference comparisons (the "naive" alternatives)
Kn_soft  = (l2+2*m2)/h2;                 % soft-layer-only
Kn_stiff = (l1+2*m1)/h1;                 % stiff-layer-only
[lA,mA]=lame((P.E_sei+P.E_sei2)/2,(P.nu_sei+P.nu_sei2)/2); Kn_avg=(lA+2*mA)/ep;  % naive average
Reff_th = h1/P.D_sei + h2/P.D_sei2;      % series diffusion resistance

opts = struct('Lx',12*ep,'Le',max(20*ep,6e-7),'nx',24,'nye',10,'nys',8,'nyel',10,'h1_frac',h1f);
opts.Kn = Kn_th; opts.Kt = Kt_th; opts.beta_eff = 0;
Le = opts.Le; Lx = opts.Lx;
mf = build_trilayer_mesh_2d(P,opts);   % FULL: regions 2 (inner) + 4 (outer)

% ---- normal (rollers + affine uy) : extract K_n ----
eps0=1e-3;
BCf.ux=struct('nodes',[mf.left,mf.right],'val',0);
BCf.uy=struct('nodes',[mf.bottom,mf.top],'val',[zeros(1,numel(mf.bottom)),eps0*mf.nodes(mf.top,2)']);
sF=solve_mech_2d(P,'full',BCf,opts);
Rn=sF.K*sF.uvec; sig_yy=sum(Rn(2*mf.top))/Lx;
Dn = samp(sF,mf,Lx/2,Le+ep,2)-samp(sF,mf,Lx/2,Le,2);
Kn_meas = sig_yy/Dn;

% ---- shear (periodic x) : extract K_t ----
gamma=1e-3; optsP=opts; optsP.periodic_x=true;
Htop=max(mf.nodes(:,2));
BCs.ux=struct('nodes',[mf.bottom,mf.top],'val',[zeros(1,numel(mf.bottom)),gamma*Htop*ones(1,numel(mf.top))]);
BCs.uy=struct('nodes',[mf.bottom,mf.top],'val',0);
sFs=solve_mech_2d(P,'full',BCs,optsP);
Rs=sFs.K*sFs.uvec; sig_xy=sum(Rs(2*mf.top-1))/Lx;
Dt = samp(sFs,mf,Lx/2,Le+ep,1)-samp(sFs,mf,Lx/2,Le,1);
Kt_meas = sig_xy/Dt;

fprintf('\n  Sublayers: inner E=%.0f GPa (h=%.0f%%), outer E=%.1f GPa (h=%.0f%%)\n',...
        P.E_sei/1e9,100*h1f,P.E_sei2/1e9,100*(1-h1f));
fprintf('\n  K_n :  measured (resolved bilayer) = %.3e   series formula = %.3e   (ratio %.4f)\n',...
        Kn_meas,Kn_th,Kn_meas/Kn_th);
fprintf('  K_t :  measured                    = %.3e   series formula = %.3e   (ratio %.4f)\n',...
        Kt_meas,Kt_th,Kt_meas/Kt_th);
fprintf('\n  NON-ADDITIVE insight (mechanical stiffness is controlled by the SOFT layer):\n');
fprintf('     K_n (series)        = %.3e Pa/m\n', Kn_th);
fprintf('     soft-layer only     = %.3e Pa/m   -> K_n / soft = %.3f  (close to 1)\n', Kn_soft, Kn_th/Kn_soft);
fprintf('     stiff-layer only    = %.3e Pa/m   -> K_n / stiff = %.4f (tiny)\n', Kn_stiff, Kn_th/Kn_stiff);
fprintf('     naive thickness-avg = %.3e Pa/m   -> K_n / avg  = %.3f  (NOT the average!)\n', Kn_avg, Kn_th/Kn_avg);
fprintf('\n  Complementary (diffusion, controlled by the BLOCKING inner layer):\n');
fprintf('     R_eff (series) = %.3e s/m ;  inner contributes %.0f%% , outer %.0f%%\n',...
        Reff_th, 100*(h1/P.D_sei)/Reff_th, 100*(h2/P.D_sei2)/Reff_th);
fprintf('\nDONE.\n');

out=struct('Kn_meas',Kn_meas,'Kn_th',Kn_th,'Kt_meas',Kt_meas,'Kt_th',Kt_th,...
   'Kn_soft',Kn_soft,'Kn_stiff',Kn_stiff,'Kn_avg',Kn_avg,'Reff_th',Reff_th,...
   'inner_R_frac',(h1/P.D_sei)/Reff_th);
end

function v=samp(sol,m,xq,yq,comp)
d=(m.nodes(:,1)-xq).^2+(m.nodes(:,2)-yq).^2; [~,k]=min(d); v=sol.u(k,comp);
end
