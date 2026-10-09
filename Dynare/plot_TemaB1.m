function fig = plot_TemaB1(experiment, p, ss, figname, pref)
% PLOT_TEMAB1  Draws the figures of Matlab/TemaB1_experimentos.m.
%
%   fig = plot_TemaB1(experiment, p, ss)
%   fig = plot_TemaB1(experiment, p, ss, figname)
%   fig = plot_TemaB1(experiment, p, ss, figname, pref)
%
%   experiment : 1,...,6, number of the experiment (and figure) in
%                TemaB1_experimentos.m
%                1 productivity shock, unanticipated      (chi = 0)
%                2 productivity shock, anticipated        (chi = 0)
%                3 productivity shock, unanticipated      (chi = 0.001)
%                4 world interest rate shock, unanticipated (chi = 0)
%                5 world interest rate shock, anticipated   (chi = 0)
%                6 world interest rate shock, unanticipated (chi = 0.001)
%   p, ss      : paths and initial steady state (see get_TemaB1_paths.m)
%   figname    : optional figure name
%   pref       : optional second set of paths (same field names), e.g. the
%                solution of TemaB1_experimentos.m, drawn as black dots on
%                top of p for comparison.
%
% Same panels, units, colors and line widths as TemaB1_experimentos.m: the
% first 51 periods are plotted (labelled 0,...,50), the red circle is the
% initial steady state and the red dashed line its level. All variables are
% in levels; as in TemaB1_experimentos.m, the interest rates are in percent
% in figures 3 and 6 and in decimal units in figures 4 and 5. Also as in
% TemaB1_experimentos.m, assets are the stock at the beginning of the period
% in figures 1, 2, 4 and 5 (fsolve) and at the end of the period in figures
% 3 and 6 (variable S of its Dynare model temaB1.mod). Titles are
% written without accents (ASCII). ylim('padded') is used where available
% (MATLAB R2021a or later), otherwise an equivalent padding is computed.

if nargin < 4 || isempty(figname)
    figname = sprintf('Dynare: Tema B1.%d', experiment);
end
if nargin < 5
    pref = [];
end

% Panels: {path field, steady state field, title, scale}
th  = {'tht',  'th0',   'Productividad (\theta_t)  [IMPULSO]',   1};
cc  = {'ct',   'css',   'Consumo (C_t)',                         1};
ll  = {'lt',   'lss',   'Trabajo (L_t)',                         1};
yy  = {'yt',   'yss',   'Producto (Y_t)',                        1};
ca  = {'cat',  'cass',  'Cuenta Corriente (CA_t)',               1};
sa  = {'st',   's0',    'Activos Externos Netos (S^*_t)',        1};
se  = {'st_end', 's0',  'Activos Externos Netos (S^*_t)',        1};
re  = {'rste', 'rlong', 'Tasa de Interes Efectiva (r^*_t)',      100};
switch experiment
    case 1
        panels = [th; cc; ll; yy; ca; sa];
        ttl = 'Choque Tecnologico Persistente No-Anticipado (\chi = 0)';
    case 2
        panels = [th; cc; ll; yy; ca; sa];
        ttl = 'Choque Tecnologico Anticipado 5 Periodos (\chi = 0)';
    case 3
        panels = [th; cc; ll; ca; se; re];
        ttl = 'Choque Tecnologico con Tasa de Interes Endogena (\chi = 0.001)';
    case 4
        panels = [{'rst', 'rlong', 'Tasa de Interes (r^*_t)  [IMPULSO]', 1}; cc; ll; yy; ca; sa];
        ttl = 'Choque a la Tasa de Interes No-Anticipado (\chi = 0)';
    case 5
        panels = [{'rst', 'rlong', 'Tasa de Interes (r^*_t)  [IMPULSO]', 1}; cc; ll; yy; ca; sa];
        ttl = 'Choque a la Tasa de Interes Anticipado 5 Periodos (\chi = 0)';
    case 6
        panels = [{'rst', 'rlong', 'Tasa Mundial (r^w_t)  [IMPULSO]', 100}; cc; ll; ca; se; re];
        ttl = 'Choque a la Tasa Mundial con Tasa Endogena (\chi = 0.001)';
    otherwise
        error('plot_TemaB1: experiment must be 1,...,6')
end

% Style of TemaB1_experimentos.m
col_dyn = [0.00, 0.38, 0.78];
col_ss  = [0.85, 0.15, 0.15];
lw_dyn  = 2.2;
lw_ss   = 1.3;
Nplot   = 51;
t_vec   = 0:(Nplot-1);

fig = figure('Name', figname, 'Position', [60, 60, 1280, 680], 'Color', 'w');
for k = 1:6
    ax = subplot(2, 3, k);
    x  = panels{k,4}*p.(panels{k,1})(1:Nplot);
    x0 = panels{k,4}*ss.(panels{k,2});
    h1 = plot(t_vec, x, 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
    h0 = plot(0, x0, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
    h2 = plot(t_vec, x0*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
    ydata = [x(:); x0];
    if ~isempty(pref)
        xr = panels{k,4}*pref.(panels{k,1})(1:Nplot);
        hr = plot(t_vec, xr, 'k.', 'MarkerSize', 10);
        ydata = [ydata; xr(:)];
    end
    hold off
    xlim([0 50]);
    padded_ylim(ax, ydata);
    grid on;
    if k == 1
        title(panels{k,3}, 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.05 0.25 0.55]);
        set(ax, 'LineWidth', 1.0, 'Box', 'on', 'FontSize', 14, 'TickDir', 'out', ...
                'XColor', [0.05 0.20 0.45], 'YColor', [0.05 0.20 0.45], ...
                'GridAlpha', 0.20, 'Color', [0.985 0.992 1.00]);
        if isempty(pref)
            legend([h1 h2], {'Trayectoria dinamica', 'Estado estacionario'}, ...
                   'Location', 'northeast', 'FontSize', 12, 'Box', 'off');
        else
            legend([h1 hr h2], {'Dynare', 'MATLAB', 'Estado estacionario'}, ...
                   'Location', 'northeast', 'FontSize', 12, 'Box', 'off');
        end
    else
        title(panels{k,3}, 'FontSize', 18, 'FontWeight', 'bold');
        set(ax, 'LineWidth', 1.0, 'Box', 'on', 'FontSize', 14, 'TickDir', 'out', ...
                'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2], ...
                'GridAlpha', 0.20, 'Color', 'w');
    end
    if k >= 4
        xlabel('Periodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
    end
end

if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
    sgtitle({ttl, ' '}, 'FontSize', 24, 'FontWeight', 'bold');
else
    % No sgtitle (Octave, MATLAB < R2018b): leave room above the panels
    axs = findobj(fig, 'Type', 'axes', '-not', 'Tag', 'legend');
    for k = 1:numel(axs)
        pos = get(axs(k), 'Position');
        set(axs(k), 'Position', [pos(1), 0.90*pos(2), pos(3), 0.90*pos(4)]);
    end
    annotation('textbox', [0 0.93 1 0.06], 'String', ttl, 'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');
end
end


function padded_ylim(ax, ydata)
% ylim('padded') where available, otherwise a 7% margin around the data
try
    ylim(ax, 'padded');
catch
    lo = min(ydata);
    hi = max(ydata);
    pad = 0.07*(hi-lo);
    if pad <= 0
        pad = 0.07*max(abs(hi), 1);
    end
    ylim(ax, [lo-pad, hi+pad]);
end
end
