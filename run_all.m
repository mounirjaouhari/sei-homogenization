function run_all(mode)
%RUN_ALL  Master driver for the SEI chemo-mechanical homogenization code.
%
%   run_all          run the full reproducible pipeline (1D + scaled + 2D)
%   run_all('1d')    1D validation only (Section 6.2-6.5, ~45 s)
%   run_all('2d')    2D anisotropy + 2D convergence only (~30 s)
%   run_all('scaled')scaled-modulus regime study only
%
%   All MATLAB, no commercial FE. Results are written to figures/ and results/.
%
%   Reproduces the numerical validation of:
%     "Asymptotic Derivation of Chemo-Mechanical Jump Conditions for the Solid
%      Electrolyte Interphase in Lithium-Ion Batteries"
%
%   See also: tests/main_validation, tests/test_2d_anisotropy

if nargin < 1, mode = 'all'; end
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root,'src')));
addpath(fullfile(root,'tests'));

fprintf('================================================================\n');
fprintf('  SEI homogenization -- reproducible pipeline   [mode = %s]\n', mode);
fprintf('================================================================\n');
t0 = tic; S = struct();

if any(strcmp(mode,{'all','1d'}))
    fprintf('\n>>> PHASE 1 : 1-D validation (cases 1-4)\n');
    main_validation();                       % steady, transient, fracture, convergence
end

if any(strcmp(mode,{'all','scaled'}))
    fprintf('\n>>> PHASE 2 : scaled-modulus regime study\n');
    S.scaled = test_scaled_modulus();
end

if any(strcmp(mode,{'all','2d'}))
    fprintf('\n>>> PHASE 3 : 2-D plane-strain anisotropy (K_n / K_t)\n');
    S.aniso = test_2d_anisotropy();
    fprintf('\n>>> PHASE 4 : 2-D convergence (normal + shear)\n');
    S.conv2d = test_2d_convergence();
    fprintf('\n>>> PHASE 5 : 2-D bilayer interface (inorganic + organic)\n');
    S.bilayer = test_2d_bilayer();
    fprintf('\n>>> PHASE 6 : curved interface (axisymmetric, curvature eps/R)\n');
    S.curved = test_2d_curved();
end

fprintf('\n================================================================\n');
fprintf('  SUMMARY\n');
fprintf('================================================================\n');
if isfield(S,'aniso')
    a = S.aniso;
    fprintf('  2D K_n/K_t (anisotropy)  : measured %.3f  vs theory %.3f  -> %s\n', ...
        a.ratio_meas, a.ratio_th, passfail(abs(a.ratio_meas/a.ratio_th-1)<0.05));
end
if isfield(S,'conv2d')
    c = S.conv2d;
    fprintf('  2D convergence rate      : normal %.2f , shear %.2f  -> %s\n', ...
        c.rateN, c.rateT, passfail(c.rateN>0.8 && c.rateT>0.8));
end
if isfield(S,'bilayer')
    b = S.bilayer;
    fprintf('  2D bilayer K_n (series)  : measured/formula %.3f ; K_n/avg %.2f (soft-dominated) -> %s\n', ...
        b.Kn_meas/b.Kn_th, b.Kn_th/b.Kn_avg, passfail(abs(b.Kn_meas/b.Kn_th-1)<0.02));
end
if isfield(S,'curved')
    cv = S.curved;
    fprintf('  curvature correction     : err ~ (eps/R)^%.2f ; flat-limit err %.1e -> %s\n', ...
        cv.slope, cv.err(1), passfail(cv.slope>0.8 && cv.err(1)<1e-3));
end
fprintf('  Total wall-clock         : %.1f s\n', toc(t0));
fprintf('================================================================\n');
end

function s = passfail(ok), if ok, s='PASS'; else, s='FAIL'; end, end
