function fig = plot_TemaB1(experiment, p, ss, figname, pref)
% PLOT_TEMAB1  Draws the same 3x2 figure as Matlab/TemaB1.m.
%
%   fig = plot_TemaB1(experiment, p, ss)
%   fig = plot_TemaB1(experiment, p, ss, figname)
%   fig = plot_TemaB1(experiment, p, ss, figname, pref)
%
%   experiment = 1 : layout of figure(1) of TemaB1.m (exogenous interest rate)
%   experiment = 2 : layout of figure(2) of TemaB1.m (endogenous risk premium)
%   p, ss          : paths and initial steady state (see get_TemaB1_paths.m)
%   figname        : optional figure name
%   pref           : optional second set of paths (same field names), e.g.
%                    the fsolve solution of TemaB1.m, drawn as black dots on
%                    top of p for comparison.
%
% As in TemaB1.m, the first 51 periods are plotted (labelled 0,...,50), the
% red circle is the initial steady state and the red dashed line its level.
% All variables are in levels; the interest rate is in percent.

if nargin < 4 || isempty(figname)
    figname = sprintf('Dynare: experimento %d', experiment);
end
if nargin < 5
    pref = [];
end

% {path field, steady state field, axis, label, scale}
if experiment == 1
    panels = {'tht',  'th0',   [0 50 8 17],    'Productividad',    1;
              'ct',   'css',   [0 50 5 9],     'Consumo',          1;
              'lt',   'lss',   [0 50 0.6 0.8], 'Trabajo',          1;
              'yt',   'yss',   [0 50 5 12],    'Producto',         1;
              'cat',  'cass',  [0 50 -2 5],    'Cuenta Corriente', 1;
              'st',   's0',    [0 50 -10 35],  'Activos',          1};
else
    panels = {'tht',  'th0',   [0 50 8 17],    'Productividad',    1;
              'ct',   'css',   [0 50 5 9],     'Consumo',          1;
              'lt',   'lss',   [0 50 0.6 0.8], 'Trabajo',          1;
              'cat',  'cass',  [0 50 -2 5],    'Cuenta Corriente', 1;
              'st',   's0',    [0 50 -10 35],  'Activos',          1;
              'rste', 'rlong', [0 50 3 12],    'Tasa de Interes',  100};
end

t   = (1:51)-1;
fig = figure('Name', figname);
for k = 1:size(panels, 1)
    subplot(3, 2, k)
    x  = panels{k,5}*p.(panels{k,1});
    x0 = panels{k,5}*ss.(panels{k,2});
    h = plot(t, x(1:51), 'b-', 0, x0, 'ro', t, x0*ones(51,1), 'r--');
    if ~isempty(pref)
        xr = panels{k,5}*pref.(panels{k,1});
        hold on
        hr = plot(t, xr(1:51), 'k.');
        hold off
        if k == 1
            legend([h(1) hr], {'Dynare', 'MATLAB (fsolve)'})
        end
    end
    axis(panels{k,3})
    xlabel(panels{k,4})
end
