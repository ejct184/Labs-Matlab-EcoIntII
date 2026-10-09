function results = compare_matlab_dynare(save_figures)
% COMPARE_MATLAB_DYNARE  Validates the Dynare replication of the experiments
% of Matlab/TemaB1_experimentos.m.
%
%   results = compare_matlab_dynare
%   results = compare_matlab_dynare(true)   % also saves the figures as PNG
%                                           % files in the folder figures/
%
% Experiments (numbered as in TemaB1_experimentos.m):
%   1, 2, 4, 5 (chi = 0)     : SOE_B1.mod  vs fsolve on ftrans_rbc, the
%                              system of TemaB1_experimentos.m (identical to
%                              Matlab/ftransB1.m)
%   3, 6       (chi = 0.001) : SOE_B1b.mod vs the Dynare model
%                              Matlab/temaB1.mod used by TemaB1_experimentos.m
%
% For each experiment the function
%   1. runs the Dynare model (perfect foresight solver);
%   2. solves the MATLAB reference:
%      - experiments 1, 2, 4, 5: fsolve with the same paths, initial guess
%        and options as TemaB1_experimentos.m ("as run"), and again with
%        tight tolerances ("exact");
%      - experiments 3, 6: Matlab/temaB1.mod run exactly as in
%        TemaB1_experimentos.m ("as run", default Dynare tolerance) and with
%        a tight tolerance, and fsolve on a MATLAB transcription of its
%        equations (subfunction ftrans_temaB1, "exact");
%   3. checks that both programs have the same initial steady state and
%      exogenous paths, that the Dynare paths satisfy the MATLAB equilibrium
%      conditions (max|F(x_Dynare)| at machine precision), and that the IRFs
%      coincide;
%   4. draws the figures of TemaB1_experimentos.m with both solutions.
%
% Dynare (6.x; with 5.x the temaB1.mod run, which uses Dynare 6 syntax, is
% skipped) must be on the MATLAB/Octave path. Run it from the Dynare folder
% of this repository (or with that folder on the path).
%
% Note: TemaB1_experimentos.m calls fsolve with a T x 6 matrix as initial
% guess. MATLAB's fsolve works internally with x0(:), so here the systems are
% solved through a wrapper that reshapes x(:) to T x nv and returns F(:).
% This gives the same iterations in MATLAB and also runs in Octave, whose
% fsolve does not accept matrix arguments.

if nargin < 1
    save_figures = false;
end
if ~exist('dynare', 'file')
    error('Dynare is not on the path (use addpath <dynare>/matlab).')
end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'Matlab'));     % original ftransB1.m
oldfolder = cd(here);
restore = onCleanup(@() cd(oldfolder));

tol = 1e-10;           % tolerance for the replication tests

% ---------------- 1. Dynare ----------------
dynare('SOE_B1', 'noclearall');              % the driver runs in the base workspace
dyn_B1 = evalin('base', 'soe_B1');
ss_B1  = evalin('base', 'ss_B1');
dynare('SOE_B1b', 'noclearall');
dyn_B1b = evalin('base', 'soe_B1b');
ss_B1b  = evalin('base', 'ss_B1b');

% The Dynare model Matlab/temaB1.mod of TemaB1_experimentos.m
boss = [];
matlab_folder = fullfile(here, '..', 'Matlab');
if exist(fullfile(matlab_folder, 'temaB1.mod'), 'file') == 2
    boss = run_temaB1(matlab_folder);
end

titles = {'Choque tecnologico no anticipado, chi = 0 (SOE_B1.mod)', ...
          'Choque tecnologico anticipado, chi = 0 (SOE_B1.mod)', ...
          'Choque tecnologico no anticipado, chi = 0.001 (SOE_B1b.mod)', ...
          'Choque de tasa de interes no anticipado, chi = 0 (SOE_B1.mod)', ...
          'Choque de tasa de interes anticipado, chi = 0 (SOE_B1.mod)', ...
          'Choque de tasa mundial no anticipado, chi = 0.001 (SOE_B1b.mod)'};

for e = 1:6
    field = sprintf('exp%d', e);
    premium = any(e == [3 6]);
    if premium
        dyn = dyn_B1b.(field);
        ssd = ss_B1b;
    else
        dyn = dyn_B1.(field);
        ssd = ss_B1;
    end

    % ---------------- 2. MATLAB reference ----------------
    m = solve_matlab(e, dyn);
    if premium && ~isempty(boss)
        m.asrun = boss.(field);
        m.boss_tight = boss.([field '_tight']);
    end

    % ---------------- 3. Checks ----------------
    vars  = {'wt', 'lt', 'yt', 'cat', 'st', 'ct'};
    names = {'w',  'l',  'y',  'ca',  's',  'c'};
    ssf   = {'wss', 'lss', 'yss', 'cass', 's0', 'css'};
    if premium
        vars{end+1} = 'rste'; names{end+1} = 'r'; ssf{end+1} = 'rlong';
    end
    nv = numel(vars);

    ss_diff = 0;
    for k = 1:nv
        ss_diff = max(ss_diff, abs(ssd.(ssf{k}) - m.ss.(ssf{k})));
    end
    path_diff = max([abs(dyn.tht - m.tht); abs(dyn.rst - m.rst)]);

    d_exact = max_diff(dyn, m.exact, vars);
    d_asrun = NaN(nv, 1);
    if ~isempty(m.asrun)
        d_asrun = max_diff(dyn, m.asrun, vars);
    end

    pass = dyn.status == 1 && ss_diff < tol && path_diff < tol && ...
           m.F_ss < tol && m.F_dynare < tol && max(d_exact) < tol;
    if ~isempty(m.asrun)
        pass = pass && max(d_asrun) < 1e-5;     % within the reference solver tolerance
    end
    if isfield(m, 'boss_tight')
        d_boss = max_diff(dyn, m.boss_tight, vars);
        pass = pass && max(d_boss) < tol;
    end

    % ---------------- Report ----------------
    fprintf('\n%s\nExperimento %d: %s\n%s\n', repmat('=', 1, 78), e, titles{e}, repmat('=', 1, 78));
    fprintf('MATLAB reference: %s\n', m.reference);
    fprintf('Dynare solver status: %d\n', dyn.status);
    fprintf('Initial steady state, max |Dynare - MATLAB|               %9.2e\n', ss_diff);
    fprintf('Exogenous paths (theta, r*), max |Dynare - MATLAB|        %9.2e\n', path_diff);
    fprintf('Max |F| of the MATLAB system at the initial steady state  %9.2e\n', m.F_ss);
    fprintf('Max |F| of the MATLAB system at the Dynare solution       %9.2e\n', m.F_dynare);
    if ~premium
        fprintf('Max |F| of Matlab/ftransB1.m at the Dynare solution       %9.2e\n', m.F_ftransB1);
    end
    if premium && isempty(boss)
        fprintf('(Matlab/temaB1.mod not run: "as run" comparison skipped)\n');
    end
    fprintf('\nMax |Dynare - MATLAB| (levels, t = 1..%d)   as run       exact\n', numel(dyn.ct));
    for k = 1:nv
        fprintf('  %-4s                                  %9.2e   %9.2e\n', names{k}, d_asrun(k), d_exact(k));
    end
    if isfield(m, 'boss_tight')
        fprintf('temaB1.mod solved with tolf = 1e-12, max |Dynare - temaB1.mod| %9.2e\n', max(d_boss));
    end
    if pass
        fprintf('RESULT: Dynare replicates the MATLAB experiment (tolerance %g).\n', tol);
    else
        fprintf('RESULT: replication test FAILED (tolerance %g).\n', tol);
    end

    % ---------------- 4. Figures ----------------
    pref = m.asrun;
    if isempty(pref)
        pref = m.exact;
    end
    fig = plot_TemaB1(e, dyn, ssd, sprintf('Dynare vs MATLAB: Tema B1.%d', e), pref);
    if save_figures
        if ~exist(fullfile(here, 'figures'), 'dir')
            mkdir(fullfile(here, 'figures'));
        end
        print(fig, '-dpng', '-r100', fullfile(here, 'figures', sprintf('TemaB1_exp%d_vs_MATLAB.png', e)));
    end

    results.(field) = struct('dynare', dyn, 'dynare_ss', ssd, 'matlab', m, ...
        'irf_diff_asrun', d_asrun, 'irf_diff_exact', d_exact, 'pass', pass);
end
end


function d = max_diff(p, q, vars)
% Max absolute difference between two sets of paths, variable by variable
d = zeros(numel(vars), 1);
for k = 1:numel(vars)
    d(k) = max(abs(p.(vars{k}) - q.(vars{k})));
end
end


function m = solve_matlab(e, dyn)
% MATLAB reference for experiment e. Parameters, steady state and the paths
% of experiments 1, 2, 4 and 5 are copied from TemaB1_experimentos.m.

beta  = 0.95;
rlong = (1/beta) - 1;
phi   = 0.50;
th0   = 10.0;
s0    = 0.0;
chi   = 0.001;

css   = (th0/(1+phi)) + (1/(1+phi))*((1-beta)/beta)*s0;
wss   = th0;
lss   = 1 - phi*(css/wss);
yss   = th0 * lss;
cass  = 0.0;

T = 100;
opts  = optimset('Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10);   % as in TemaB1_experimentos.m
tight = optimset('Display', 'off', 'TolFun', 1e-14, 'TolX', 1e-14, ...
                 'MaxIter', 1000, 'MaxFunEvals', 1e6);

tht = th0 * ones(T, 1);
rst = rlong * ones(T, 1);
switch e
    case 1
        tht(1) = 15.0;
        for t = 1:29
            tht(t+1) = th0 + 0.88 * (tht(t) - th0);
        end
    case 2
        tht(5) = 15.0;
        for t = 5:34
            tht(t+1) = th0 + 0.88 * (tht(t) - th0);
        end
    case 4
        rst(1) = 0.10;
        for t = 1:29
            rst(t+1) = rlong + 0.88 * (rst(t) - rlong);
        end
    case 5
        rst(5) = 0.10;
        for t = 5:34
            rst(t+1) = rlong + 0.88 * (rst(t) - rlong);
        end
    case 3
        % AR(1) of temaB1.mod hit by e_theta = 5 in period 1
        tht(1) = 15.0;
        for t = 1:T-1
            tht(t+1) = th0 + 0.88 * (tht(t) - th0);
        end
    case 6
        % AR(1) of temaB1.mod hit by e_rw = 0.10 - rlong in period 1
        rst(1) = 0.10;
        for t = 1:T-1
            rst(t+1) = rlong + 0.88 * (rst(t) - rlong);
        end
end

if any(e == [1 2 4 5])
    f  = @(x) ftrans_rbc(x, T, beta, phi, s0, rst, tht);
    x0 = ones(T, 1) * [wss, lss, yss, cass, s0, css];
    m.reference = 'ftrans_rbc of TemaB1_experimentos.m (= Matlab/ftransB1.m), fsolve';
else
    f  = @(x) ftrans_temaB1(x, T, beta, phi, chi, s0, rst, tht, css, rlong);
    x0 = ones(T, 1) * [wss, lss, yss, cass, s0, css, rlong];
    m.reference = 'Matlab/temaB1.mod (as run) and its MATLAB transcription (exact)';
end
nv = size(x0, 2);
g  = @(z) reshape(f(reshape(z, T, nv)), [], 1);

m.asrun = [];
z = fsolve(g, x0(:), opts);
if any(e == [1 2 4 5])
    m.asrun = to_struct(reshape(z, T, nv), tht, rst);
end
z = fsolve(g, z, tight);
m.exact = to_struct(reshape(z, T, nv), tht, rst);

% MATLAB equilibrium conditions evaluated at the Dynare solution
ydyn = [dyn.wt dyn.lt dyn.yt dyn.cat dyn.st dyn.ct];
if nv == 7
    ydyn = [ydyn dyn.rste];
end
m.F_dynare = max(max(abs(f(ydyn))));
if nv == 6
    m.F_ftransB1 = max(max(abs(ftransB1_globals(ydyn, T, beta, phi, s0, rst, tht))));
end

% Initial steady state: solves the MATLAB system when there is no shock
if nv == 6
    F = ftrans_rbc(x0, T, beta, phi, s0, rlong*ones(T,1), th0*ones(T,1));
else
    F = ftrans_temaB1(x0, T, beta, phi, chi, s0, rlong*ones(T,1), th0*ones(T,1), css, rlong);
end
m.F_ss = max(abs(F(:)));

m.tht = tht;
m.rst = rst;
m.ss  = struct('wss', wss, 'lss', lss, 'yss', yss, 'cass', cass, ...
               's0', s0, 'css', css, 'th0', th0, 'rlong', rlong);
end


function p = to_struct(y, tht, rst)
% Columns of the fsolve solution, with the names used in this folder
p.wt  = y(:,1);
p.lt  = y(:,2);
p.yt  = y(:,3);
p.cat = y(:,4);
p.st  = y(:,5);
p.ct  = y(:,6);
if size(y, 2) == 7
    p.rste = y(:,7);
end
p.st_end = p.st + p.cat;
p.tht = tht;
p.rst = rst;
end


function F = ftrans_rbc(x, T, beta, phi, s0, rst, tht)
% Copied verbatim from the local function ftrans_rbc of
% Matlab/TemaB1_experimentos.m (identical to Matlab/ftransB1.m).
    wt  = x(:, 1);
    lt  = x(:, 2);
    yt  = x(:, 3);
    cat = x(:, 4);
    st  = x(:, 5);
    ct  = x(:, 6);

    fw  = wt - tht;
    fl  = lt - (1 - phi * (ct ./ wt));
    fy  = yt - tht .* lt;
    fca = cat - (yt - ct + rst .* st);

    fs = zeros(T, 1);
    fs(1)    = st(1) - s0;
    fs(2:T)  = st(2:T) - (st(1:T-1) + cat(1:T-1));

    fc = zeros(T, 1);
    fc(1:T-1) = ct(2:T) - beta * (1 + rst(2:T)) .* ct(1:T-1);
    fc(T)     = st(T-1) - st(T);

    F = [fw, fl, fy, fca, fs, fc];
end


function F = ftransB1_globals(x, T_, beta_, phi_, s0_, rst_, tht_)
% Evaluates the original Matlab/ftransB1.m (which reads globals)
global beta phi s0
global T rst tht
beta = beta_; phi = phi_; s0 = s0_; T = T_; rst = rst_; tht = tht_;
F = ftransB1(x);
end


function F = ftrans_temaB1(x, T, beta, phi, chi, s0, rw, tht, css, rlong)
% MATLAB transcription of the equations of temaB1.mod (the Dynare model used
% by TemaB1_experimentos.m for experiments 3 and 6), written for this
% validation with s(t) = stock of assets at the beginning of period t
% (S(t-1) in temaB1.mod):
%   w = theta, L = 1 - phi*C/w, Y = theta*L, CA = Y - C + r*s,
%   s(t+1) = s(t) + CA(t), s(1) = s0, C(t+1) = beta*(1+r(t+1))*C(t),
%   r = rw - chi*s, and Dynare's terminal condition: C and r equal to their
%   steady state values in period T+1.
    wt  = x(:, 1);
    lt  = x(:, 2);
    yt  = x(:, 3);
    cat = x(:, 4);
    st  = x(:, 5);
    ct  = x(:, 6);
    rt  = x(:, 7);

    fw  = wt - tht;
    fl  = lt - (1 - phi * (ct ./ wt));
    fy  = yt - tht .* lt;
    fca = cat - (yt - ct + rt .* st);

    fs = zeros(T, 1);
    fs(1)    = st(1) - s0;
    fs(2:T)  = st(2:T) - (st(1:T-1) + cat(1:T-1));

    fc = zeros(T, 1);
    fc(1:T-1) = ct(2:T) - beta * (1 + rt(2:T)) .* ct(1:T-1);
    fc(T)     = css - beta * (1 + rlong) * ct(T);

    fr = rt - (rw - chi * st);

    F = [fw, fl, fy, fca, fs, fc, fr];
end


function out = run_temaB1(matlab_folder)
% Runs experiments 3 and 6 with Matlab/temaB1.mod exactly as
% TemaB1_experimentos.m does (Dynare 6 syntax), with the default tolerance
% and with tolf = 1e-12. Dynare needs the .mod file in the current folder.
out = [];
oldfolder = cd(matlab_folder);
restore = onCleanup(@() cd(oldfolder));
try
    evalc('dynare temaB1.mod noclearall');
    base_cmd = ['options_.periods = 100; set_param_value(''chi'', 0.001); ' ...
                'oo_ = perfect_foresight_setup(M_, options_, oo_); ' ...
                '[oo_, ~] = perfect_foresight_solver(M_, options_, oo_);'];
    shock3 = 'M_.det_shocks = struct(''exo_det'', false, ''exo_id'', 1, ''type'', ''level'', ''periods'', 1:1, ''value'', 5.0); ';
    shock6 = 'M_.det_shocks = struct(''exo_det'', false, ''exo_id'', 2, ''type'', ''level'', ''periods'', 1:1, ''value'', 0.10 - ((1/0.95) - 1)); ';
    tolcmd = {'', 'options_.dynatol.f = 1e-12; options_.dynatol.x = 1e-12; '};
    suffix = {'', '_tight'};
    for k = 1:2
        evalin('base', [tolcmd{k} shock3 base_cmd]);
        out.(['exp3' suffix{k}]) = temaB1_paths();
        evalin('base', [tolcmd{k} shock6 base_cmd]);
        out.(['exp6' suffix{k}]) = temaB1_paths();
    end
catch err
    fprintf('\nCould not run temaB1.mod (%s); "as run" comparison skipped.\n', err.message);
    out = [];
end
end


function p = temaB1_paths()
% Paths of temaB1.mod (periods 1,...,T), with the names used in this folder
M  = evalin('base', 'M_');
oo = evalin('base', 'oo_');
X  = oo.endo_simul;
T  = size(X, 2) - 2;
v  = @(name) X(strcmp(name, cellstr(M.endo_names)), :)';
S  = v('S');
p.wt     = v('w');     p.wt  = p.wt(2:T+1);
p.lt     = v('L');     p.lt  = p.lt(2:T+1);
p.yt     = v('Y');     p.yt  = p.yt(2:T+1);
p.cat    = v('CA');    p.cat = p.cat(2:T+1);
p.ct     = v('C');     p.ct  = p.ct(2:T+1);
p.rste   = v('r');     p.rste = p.rste(2:T+1);
p.st     = S(1:T);     % beginning of period stock = S(t-1)
p.st_end = S(2:T+1);   % end of period stock S(t)
p.tht    = v('theta'); p.tht = p.tht(2:T+1);
p.rst    = v('rw');    p.rst = p.rst(2:T+1);
end
