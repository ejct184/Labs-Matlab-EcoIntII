function [p, ss] = get_TemaB1_paths(M_, oo_)
% GET_TEMAB1_PATHS  Extracts a Dynare perfect foresight solution in the format
% used by Matlab/TemaB1.m.
%
%   [p, ss] = get_TemaB1_paths(M_, oo_)
%
% Inputs : M_, oo_  Dynare structures after perfect_foresight_solver.
% Outputs: p   struct with the T x 1 transition paths (periods 1,...,T, i.e.
%              the rows of the fsolve solution in TemaB1.m):
%                p.wt, p.lt, p.yt, p.cat, p.st, p.ct  (both models)
%                p.rste  (endogenous interest rate, model with risk premium)
%                p.rst   (world interest rate path)
%                p.tht   (productivity path)
%          ss  struct with the initial steady state, named as in TemaB1.m:
%                ss.wss, ss.lss, ss.yss, ss.cass, ss.s0, ss.css, ss.th0,
%                ss.rlong

nlag  = M_.maximum_lag;
T     = size(oo_.endo_simul, 2) - nlag - M_.maximum_lead;
rows  = nlag + (1:T);                 % columns of periods 1,...,T

en    = cellstr(M_.endo_names);       % cell array in Dynare >= 4.6
xn    = cellstr(M_.exo_names);
endo  = @(name) oo_.endo_simul(strcmp(name, en), rows)';
exo   = @(name) oo_.exo_simul(rows, strcmp(name, xn));
endss = @(name) oo_.steady_state(strcmp(name, en));
exoss = @(name) oo_.exo_steady_state(strcmp(name, xn));

p.wt  = endo('w');
p.lt  = endo('l');
p.yt  = endo('y');
p.cat = endo('ca');
p.st  = endo('s');
p.ct  = endo('c');
p.tht = exo('theta');
p.rst = exo('rst');
if any(strcmp('r', en))
    p.rste = endo('r');
end

ss.wss   = endss('w');
ss.lss   = endss('l');
ss.yss   = endss('y');
ss.cass  = endss('ca');
ss.s0    = endss('s');
ss.css   = endss('c');
ss.th0   = exoss('theta');
ss.rlong = exoss('rst');
