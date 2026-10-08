function results = compare_matlab_dynare(save_figures)
% COMPARE_MATLAB_DYNARE  Validates the Dynare replication of Matlab/TemaB1.m.
%
%   results = compare_matlab_dynare
%   results = compare_matlab_dynare(true)   % also saves the figures as PNG
%                                           % files in the folder figures/
%
% For each experiment of TemaB1.m (I: exogenous interest rate, SOE_B1.mod;
% II: endogenous risk premium, SOE_B1b.mod) the function
%   1. runs the Dynare model (perfect foresight solver);
%   2. solves the ORIGINAL MATLAB system, calling Matlab/ftransB1.m and
%      Matlab/ftransB1b.m unchanged, with the same globals, initial guess and
%      fsolve options as TemaB1.m ("as run"), and again with tight
%      tolerances (exact solution of the same system);
%   3. checks that
%        a. both programs have the same initial steady state, and that it
%           solves the original MATLAB system when there is no shock;
%        b. the Dynare paths satisfy the original MATLAB equilibrium
%           conditions, i.e. max|F(x_Dynare)| is at machine precision;
%        c. the IRFs (levels, t = 1,...,T) coincide;
%   4. draws the figures of TemaB1.m with both solutions superimposed.
%
% Dynare (5.x or 6.x) must be on the MATLAB/Octave path. Run it from the
% Dynare folder of this repository (or with that folder on the path).
%
% Note: the original call fsolve(@ftransB1, x0) uses a T x 6 matrix x0.
% MATLAB's fsolve works internally with x0(:), so here the original
% functions are called through a wrapper that reshapes x(:) to T x nv and
% returns F(:). This gives the same iterations in MATLAB and also runs in
% Octave, whose fsolve does not accept matrix arguments.

if nargin < 1
    save_figures = false;
end
if ~exist('dynare', 'file')
    error('Dynare is not on the path (use addpath <dynare>/matlab).')
end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'Matlab'));     % original ftransB1.m, ftransB1b.m
oldfolder = cd(here);
restore = onCleanup(@() cd(oldfolder));

tol = 1e-10;           % tolerance for the replication tests

mods   = {'SOE_B1', 'SOE_B1b'};
titles = {'Experimento I: tasa de interes exogena (ftransB1.m vs SOE_B1.mod)', ...
          'Experimento II: prima de riesgo endogena (ftransB1b.m vs SOE_B1b.mod)'};

for e = 1:2
    % ---------------- 1. Dynare ----------------
    dynare(mods{e}, 'noclearall');          % the driver runs in the base workspace
    M  = evalin('base', 'M_');
    oo = evalin('base', 'oo_');
    [dyn, ssd] = get_TemaB1_paths(M, oo);
    dyn_status = oo.deterministic_simulation.status;

    % ---------------- 2. Original MATLAB code ----------------
    m = solve_TemaB1(e, dyn);

    % ---------------- 3. Checks ----------------
    vars  = {'wt', 'lt', 'yt', 'cat', 'st', 'ct'};
    names = {'w',  'l',  'y',  'ca',  's',  'c'};
    ssf   = {'wss', 'lss', 'yss', 'cass', 's0', 'css'};
    if e == 2
        vars{end+1} = 'rste'; names{end+1} = 'r'; ssf{end+1} = 'rlong';
    end
    nv = numel(vars);

    ss_dyn = zeros(nv, 1);
    ss_mat = zeros(nv, 1);
    for k = 1:nv
        ss_dyn(k) = ssd.(ssf{k});
        ss_mat(k) = m.ss.(ssf{k});
    end

    T = numel(dyn.ct);
    d_asrun = zeros(nv, 2);       % columns: t = 1..T, t = 1..51 (plotted)
    d_tight = zeros(nv, 2);
    for k = 1:nv
        a = abs(dyn.(vars{k}) - m.asrun.(vars{k}));
        b = abs(dyn.(vars{k}) - m.tight.(vars{k}));
        d_asrun(k,:) = [max(a) max(a(1:51))];
        d_tight(k,:) = [max(b) max(b(1:51))];
    end
    d_tht = max(abs(dyn.tht - m.tht));

    pass = dyn_status == 1 && max(abs(ss_dyn - ss_mat)) < tol && ...
           m.F_ss < tol && m.F_dynare < tol && max(d_tight(:,1)) < tol && ...
           d_tht < tol;

    % ---------------- Report ----------------
    fprintf('\n%s\n%s\n%s\n', repmat('=', 1, 78), titles{e}, repmat('=', 1, 78));
    fprintf('Dynare solver status: %d    fsolve exitflag (TemaB1.m options): %d    (tight): %d\n', ...
            dyn_status, m.exitflag_asrun, m.exitflag_tight);
    fprintf('\nInitial steady state     TemaB1.m            Dynare              |diff|\n');
    for k = 1:nv
        fprintf('  %-4s              %18.15f  %18.15f  %9.2e\n', names{k}, ss_mat(k), ss_dyn(k), abs(ss_dyn(k) - ss_mat(k)));
    end
    fprintf('\nProductivity path: max|theta_Dynare - tht_TemaB1| = %9.2e\n', d_tht);
    fprintf('\nMax |F| of the original MATLAB system (%s) evaluated at:\n', func2str(m.f));
    fprintf('  initial steady state, no shock         %9.2e\n', m.F_ss);
    fprintf('  Dynare solution                        %9.2e\n', m.F_dynare);
    fprintf('  fsolve solution, TemaB1.m options      %9.2e\n', m.F_asrun);
    fprintf('  fsolve solution, tight tolerances      %9.2e\n', m.F_tight);
    fprintf('\nMax |Dynare - MATLAB| (levels)  vs fsolve (TemaB1.m options) |  vs fsolve (tight tol.)\n');
    fprintf('                     t = 1..%d     t = 1..51 |  t = 1..%d     t = 1..51\n', T, T);
    for k = 1:nv
        fprintf('  %-4s              %9.2e    %9.2e  |  %9.2e    %9.2e\n', names{k}, d_asrun(k,1), d_asrun(k,2), d_tight(k,1), d_tight(k,2));
    end
    if pass
        fprintf('\nRESULT: Dynare replicates the MATLAB model (tolerance %g).\n', tol);
    else
        fprintf('\nRESULT: replication test FAILED (tolerance %g).\n', tol);
    end

    % ---------------- 4. Figures ----------------
    fig = plot_TemaB1(e, dyn, ssd, sprintf('Dynare vs MATLAB: experimento %d', e), m.asrun);
    if save_figures
        if ~exist(fullfile(here, 'figures'), 'dir')
            mkdir(fullfile(here, 'figures'));
        end
        print(fig, '-dpng', '-r100', fullfile(here, 'figures', sprintf('%s_vs_MATLAB.png', mods{e})));
    end

    results.(mods{e}) = struct('dynare', dyn, 'dynare_ss', ssd, 'matlab', m, ...
        'ss_diff', abs(ss_dyn - ss_mat), 'irf_diff_asrun', d_asrun, ...
        'irf_diff_tight', d_tight, 'pass', pass);
end
end


function m = solve_TemaB1(experiment, dyn)
% Solves experiment 1 or 2 of TemaB1.m with the original ftransB1.m or
% ftransB1b.m. The setup below is copied from TemaB1.m.

global beta phi s0 chi
global T rst tht

beta  = 0.95;
rlong = (1/beta)-1;
phi   = 0.5;
s0    = 0;
th0   = 10;
chi   = 0.001;

css   = (th0/(1+phi))+(1/(1+phi))*((1-beta)/beta)*s0;
wss   = th0;
lss   = 1-phi*(css/wss);
yss   = th0*lss;
cass  = 0;

T   = 100;
rst = rlong*ones(T,1);
tht = th0*ones(T,1);
tht(1) = 15;
for i = 1:29
    tht(i+1) = th0+0.88*(tht(i)-th0);
end

if experiment == 1
    f  = @ftransB1;
    ss = [wss; lss; yss; cass; s0; css];
else
    f  = @ftransB1b;
    ss = [wss; lss; yss; cass; s0; css; rlong];
end
nv = numel(ss);
x0 = ones(1,T)'*ss';
g  = @(z) reshape(f(reshape(z, T, nv)), [], 1);

options = optimset('Display', 'off');                    % as in TemaB1.m
[z, Fv, m.exitflag_asrun] = fsolve(g, x0(:), options);
m.F_asrun = max(abs(Fv));
yasrun = reshape(z, T, nv);

options = optimset('Display', 'off', 'TolFun', 1e-14, 'TolX', 1e-14, ...
                   'MaxIter', 1000, 'MaxFunEvals', 1e6);
[z, Fv, m.exitflag_tight] = fsolve(g, x0(:), options);
m.F_tight = max(abs(Fv));
ytight = reshape(z, T, nv);

% Original equilibrium conditions evaluated at the Dynare solution
ydyn = [dyn.wt dyn.lt dyn.yt dyn.cat dyn.st dyn.ct];
if experiment == 2
    ydyn = [ydyn dyn.rste];
end
m.F_dynare = max(max(abs(f(ydyn))));

% Initial steady state: solves the original system when there is no shock
tht_shock = tht;
tht = th0*ones(T,1);
m.F_ss = max(max(abs(f(x0))));
tht = tht_shock;

m.f     = f;
m.tht   = tht;
m.asrun = to_struct(yasrun, tht);
m.tight = to_struct(ytight, tht);
m.ss    = struct('wss', wss, 'lss', lss, 'yss', yss, 'cass', cass, ...
                 's0', s0, 'css', css, 'th0', th0, 'rlong', rlong);
end


function p = to_struct(y, tht)
% Columns of the fsolve solution, named as in TemaB1.m
p.wt  = y(:,1);
p.lt  = y(:,2);
p.yt  = y(:,3);
p.cat = y(:,4);
p.st  = y(:,5);
p.ct  = y(:,6);
if size(y, 2) == 7
    p.rste = y(:,7);
end
p.tht = tht;
end
