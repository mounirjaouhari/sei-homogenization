function make_profiles_2d()
%MAKE_PROFILES_2D  Solution-profile snapshots: resolved (FULL) vs homogenized (EFF).
%
%   Produces, for a plane-strain trilayer loaded in (a) uniaxial normal strain and
%   (b) simple shear:
%     - figures/fig_profiles_2d.png : through-thickness line profiles u(y) at mid-span,
%       FULL (resolved, SEI explicitly meshed) overlaid with EFF (zero-thickness
%       interface). The displacement JUMP across the SEI / interface is visible.
%     - figures/fig_fields_2d_normal.png, fig_fields_2d_shear.png : 2-D color maps
%       of the relevant displacement component on the (amplified) deformed mesh,
%       FULL vs EFF side by side.
%
%   All from the MATLAB 2-D solver (src/solvers/solve_mech_2d.m).

addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src')));
figdir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'figures');
if ~exist(figdir,'dir'), mkdir(figdir); end

P = config_parameters('intermediate');
P.E2=P.E1; P.nu2=P.nu1; P.beta2=0;                 % symmetric bulk, isolate SEI
lam=P.E_sei*P.nu_sei/((1+P.nu_sei)*(1-2*P.nu_sei)); mu=P.E_sei/(2*(1+P.nu_sei));

Le = 6e-7;  P.eps = 0.12*Le;  ep=P.eps;            % thick-ish SEI: visible in the plot
opts = struct('Lx',6*ep,'Le',Le,'nx',16,'nye',14,'nys',8,'nyel',14);
opts.Kn=(lam+2*mu)/ep; opts.Kt=mu/ep; opts.beta_eff=0;
Lx=opts.Lx;
mf=build_trilayer_mesh_2d(P,opts); me=build_eff_mesh_2d(P,opts);

% ----- solve : normal (rollers + affine uy) and shear (periodic-x) -----
eps0=2e-2; gamma=2e-2;
bcN=@(m) struct('ux',struct('nodes',[m.left,m.right],'val',0), ...
   'uy',struct('nodes',[m.bottom,m.top],'val',[zeros(1,numel(m.bottom)),eps0*m.nodes(m.top,2)']));
sFn=solve_mech_2d(P,'full',bcN(mf),opts);  sEn=solve_mech_2d(P,'eff',bcN(me),opts);
optsP=opts; optsP.periodic_x=true;
bcS=@(m) struct('ux',struct('nodes',[m.bottom,m.top], ...
   'val',[zeros(1,numel(m.bottom)),gamma*max(m.nodes(:,2))*ones(1,numel(m.top))]), ...
   'uy',struct('nodes',[m.bottom,m.top],'val',0));
sFs=solve_mech_2d(P,'full',bcS(mf),optsP); sEs=solve_mech_2d(P,'eff',bcS(me),optsP);

% ================= line profiles u(y) at mid-span =================
[yF_n,uF_n]=profile(sFn,mf,Lx/2,2);  [yE_n,uE_n]=profile(sEn,me,Lx/2,2);
[yF_s,uF_s]=profile(sFs,mf,Lx/2,1);  [yE_s,uE_s]=profile(sEs,me,Lx/2,1);

f=figure('Position',[100 100 950 420],'Color','w','MenuBar','none','ToolBar','none');
ax1=subplot(1,2,1); oneprofile(ax1,uF_n,yF_n,uE_n,yE_n,'u_y (nm)','(a) normal mode',Le,ep);
ax2=subplot(1,2,2); oneprofile(ax2,uF_s,yF_s,uE_s,yE_s,'u_x (nm)','(b) shear mode',Le,ep);
exportgraphics(f, fullfile(figdir,'fig_profiles_2d.png'),'Resolution',200); close(f);
fprintf('  wrote fig_profiles_2d.png\n');

% ================= 2-D field maps (deformed) =================
fieldmap(sFn,sEn,2,80, 'u_y  (normal mode)', fullfile(figdir,'fig_fields_2d_normal.png'),Le,ep);
fieldmap(sFs,sEs,1,80, 'u_x  (shear mode)',  fullfile(figdir,'fig_fields_2d_shear.png'), Le,ep);
fprintf('  wrote fig_fields_2d_normal.png, fig_fields_2d_shear.png\n');
fprintf('DONE.\n');
end

% ---------- helpers ----------
function [y,u]=profile(sol,m,xq,comp)
tol=0.5*(m.xv(2)-m.xv(1));
idx=find(abs(m.nodes(:,1)-xq)<tol);
[y,o]=sort(m.nodes(idx,2)); u=sol.u(idx(o),comp);
end
function oneprofile(ax, uF,yF, uE,yE, xlab, ttl, Le,ep)
hold(ax,'on'); box(ax,'on'); ax.Toolbar.Visible='off';
pF=plot(ax, uF*1e9, yF*1e9,'b-','LineWidth',2);
pE=plot(ax, uE*1e9, yE*1e9,'r--','LineWidth',2);
xl=xlim(ax); yb=[Le Le+ep]*1e9;
pb=patch(ax,[xl(1) xl(2) xl(2) xl(1)],[yb(1) yb(1) yb(2) yb(2)],[0.8 0.8 0.85], ...
   'EdgeColor','none','FaceAlpha',0.55); uistack(pb,'bottom');
xlabel(ax,xlab); ylabel(ax,'through-thickness  y (nm)'); title(ax,ttl);
legend(ax,[pF pE pb],{'FULL (resolved)','EFF (homogenized)','SEI band'}, ...
   'Location','southeast','FontSize',8);
end
function fieldmap(solF,solE,comp,scale,ttl,fname,Le,ep)
f=figure('Position',[100 100 900 430],'Color','w');
ax1=subplot(1,2,1); drawfield(solF,comp,scale,ax1); title(sprintf('FULL (resolved):  %s',ttl));
yline(Le*1e9,'k:'); yline((Le+ep)*1e9,'k:');
ax2=subplot(1,2,2); drawfield(solE,comp,scale,ax2); title(sprintf('EFF (homogenized):  %s',ttl));
yline(Le*1e9,'k:');
exportgraphics(f,fname,'Resolution',200); close(f);
end
function drawfield(sol,comp,scale,ax)
m=sol.m; X=(m.nodes(:,1)+scale*sol.u(:,1))*1e9; Y=(m.nodes(:,2)+scale*sol.u(:,2))*1e9;
C=sol.u(:,comp)*1e9;
patch('Faces',m.elems,'Vertices',[X Y],'FaceVertexCData',C, ...
   'FaceColor','interp','EdgeColor',[0.6 0.6 0.6],'EdgeAlpha',0.25,'Parent',ax);
axis(ax,'equal','tight'); colormap(ax,parula); cb=colorbar(ax); cb.Label.String='nm';
xlabel(ax,'x (nm)'); ylabel(ax,'y (nm)');
end
