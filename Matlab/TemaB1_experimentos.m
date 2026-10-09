% =========================================================================
% ECONOMIA INTERNACIONAL II - SIMULACIONES TEMA B.1
% Profesor: Carlos Urrutia
%
% Este script reproduce todos los graficos del Tema B.1 de las notas de clase:
%   - Figura 1 : Choque tecnologico persistente no-anticipado (chi = 0)
%   - Figura 2 : Choque tecnologico anticipado (se anuncia en t=1, impacta en t=5)
%   - Figura 3 : Choque tecnologico no-anticipado con tasa endogena (chi = 0.001)
%   - Figura 4 : Choque de tasa de interes persistente no-anticipado
%   - Figura 5 : Choque de tasa de interes anticipado
%   - Figura 6 : Choque de tasa de interes con tasa endogena (chi = 0.001)
%
% =========================================================================

clear; close all; clc;

fprintf('=======================================================\n');
fprintf('  Iniciando simulaciones numericas del Tema B.1\n');
fprintf('=======================================================\n\n');

% -------------------------------------------------------------------------
% 1. Parametros y Estado Estacionario Inicial
% -------------------------------------------------------------------------
beta  = 0.95;           % Factor de descuento intertemporal
rlong = (1/beta) - 1;   % Tasa de interes real de largo plazo (~5.263%)
phi   = 0.50;           % Peso del ocio en la funcion de utilidad
th0   = 10.0;           % Productividad inicial
s0    = 0.0;            % Activos internacionales netos iniciales
chi   = 0.001;          % Sensibilidad de la prima de riesgo por deuda

% Estado estacionario analitico (diapositiva 7):
css   = (th0/(1+phi)) + (1/(1+phi))*((1-beta)/beta)*s0;  % Consumo = 6.6667
wss   = th0;                                            % Salario = 10
lss   = 1 - phi*(css/wss);                              % Trabajo = 0.6667
yss   = th0 * lss;                                      % Producto = 6.6667
cass  = 0.0;                                            % Cuenta corriente = 0

fprintf('Estado estacionario inicial:\n');
fprintf('  theta = %.2f,  r* = %.4f (%.2f%%)\n', th0, rlong, rlong*100);
fprintf('  w     = %.2f,  C  = %.4f,  L = %.4f,  Y = %.4f,  CA = %.2f,  S = %.2f\n\n', ...
        wss, css, lss, yss, cass, s0);

% Configuracion de simulacion y graficacion
T     = 100;    % Periodos de simulacion
Nplot = 51;     % Periodos a graficar: t = 0, 1, ..., 50
t_vec = 0:(Nplot-1);

% Estilos de diseno
col_dyn = [0.00, 0.38, 0.78];   % Azul cobalto para trayectoria dinamica
col_ss  = [0.85, 0.15, 0.15];   % Rojo carmesi para estado estacionario
lw_dyn  = 2.2;                  % Grosor de linea de transicion
lw_ss   = 1.3;                  % Grosor de linea de estado estacionario
fig_pos = [60, 60, 1280, 680];  % Dimension 16:9 ideal para diapositivas landscape

opts = optimset('Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10);

% Compilar el modelo Dynare base
fprintf('Compilando modelo en Dynare (temaB1.mod)...\n');
evalc('dynare temaB1.mod noclearall');
fprintf('Modelo Dynare compilado exitosamente.\n\n');

%% =========================================================================
% EXPERIMENTO 1: Choque tecnologico persistente y no-anticipado
% =========================================================================
fprintf('Ejecutando Experimento 1: Choque tecnologico no-anticipado...\n');

tht_exp1 = th0 * ones(T, 1);
tht_exp1(1) = 15.0;
for t = 1:29
    tht_exp1(t+1) = th0 + 0.88 * (tht_exp1(t) - th0);
end
rst_exp1 = rlong * ones(T, 1);

x0_exp1 = ones(T, 1) * [wss, lss, yss, cass, s0, css];
sol_exp1 = fsolve(@(x) ftrans_rbc(x, T, beta, phi, s0, rst_exp1, tht_exp1), x0_exp1, opts);

wt_1  = sol_exp1(:, 1);
lt_1  = sol_exp1(:, 2);
yt_1  = sol_exp1(:, 3);
cat_1 = sol_exp1(:, 4);
st_1  = sol_exp1(:, 5);
ct_1  = sol_exp1(:, 6);

fig1 = figure('Name', 'Tema B1.1', 'Position', fig_pos, 'Color', 'w');

% Fila 1: Impulso (Productividad), Consumo, Trabajo
ax1 = subplot(2, 3, 1);
plot(t_vec, tht_exp1(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(t_vec, th0*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
plot(0, th0, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
xlim([0 50]);ylim("padded");grid on;
title('Productividad (\theta_{t})  [IMPULSO]', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.05 0.25 0.55]);
formatear_eje_impulso(ax1);
legend({'Trayectoria dinámica', 'Estado estacionario'}, 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');

ax2 = subplot(2, 3, 2);
plot(t_vec, ct_1(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, css, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, css*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Consumo (C_t)', 'FontSize', 18, 'FontWeight', 'bold');
formatear_eje(ax2);

ax3 = subplot(2, 3, 3);
plot(t_vec, lt_1(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, lss, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, lss*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Trabajo (L_t)', 'FontSize', 18, 'FontWeight', 'bold');
formatear_eje(ax3);

% Fila 2: Producto, Cuenta Corriente, Activos
ax4 = subplot(2, 3, 4);
plot(t_vec, yt_1(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, yss, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, yss*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Producto (Y_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax4);

ax5 = subplot(2, 3, 5);
plot(t_vec, cat_1(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, cass, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, cass*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Cuenta Corriente (CA_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax5);

ax6 = subplot(2, 3, 6);
plot(t_vec, st_1(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, s0, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, s0*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Activos Externos Netos (S^*_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax6);

sgtitle({'Choque Tecnológico Persistente No-Anticipado (\chi = 0)',' '}, 'FontSize', 24, 'FontWeight', 'bold');
exportgraphics(fig1, 'FigB1.1.pdf')
fprintf('  -> Guardada: Figura_B1_Exp1_Tecnologico_NoAnticipado.png\n\n');


%% =========================================================================
% EXPERIMENTO 2: Choque tecnologico anticipado 
% =========================================================================
fprintf('Ejecutando Experimento 2: Choque tecnologico anticipado ...\n');

tht_exp2 = th0 * ones(T, 1);
tht_exp2(5) = 15.0;
for t = 5:34
    tht_exp2(t+1) = th0 + 0.88 * (tht_exp2(t) - th0);
end
rst_exp2 = rlong * ones(T, 1);

x0_exp2 = ones(T, 1) * [wss, lss, yss, cass, s0, css];
sol_exp2 = fsolve(@(x) ftrans_rbc(x, T, beta, phi, s0, rst_exp2, tht_exp2), x0_exp2, opts);

wt_2  = sol_exp2(:, 1);
lt_2  = sol_exp2(:, 2);
yt_2  = sol_exp2(:, 3);
cat_2 = sol_exp2(:, 4);
st_2  = sol_exp2(:, 5);
ct_2  = sol_exp2(:, 6);

fig2 = figure('Name', 'Tema B1.2', 'Position', fig_pos, 'Color', 'w');

ax1 = subplot(2, 3, 1);
plot(t_vec, tht_exp2(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(t_vec, th0*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
plot(0, th0, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
xlim([0 50]);ylim("padded"); grid on;
title('Productividad (\theta_t)  [IMPULSO]', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.05 0.25 0.55]);
formatear_eje_impulso(ax1);
legend({'Trayectoria dinámica', 'Estado estacionario'}, 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');

ax2 = subplot(2, 3, 2);
plot(t_vec, ct_2(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, css, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, css*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Consumo (C_t)', 'FontSize', 18, 'FontWeight', 'bold');
formatear_eje(ax2);

ax3 = subplot(2, 3, 3);
plot(t_vec, lt_2(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, lss, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, lss*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Trabajo (L_t)', 'FontSize', 18, 'FontWeight', 'bold');
formatear_eje(ax3);

ax4 = subplot(2, 3, 4);
plot(t_vec, yt_2(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, yss, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, yss*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Producto (Y_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax4);

ax5 = subplot(2, 3, 5);
plot(t_vec, cat_2(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, cass, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, cass*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Cuenta Corriente (CA_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax5);

ax6 = subplot(2, 3, 6);
plot(t_vec, st_2(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, s0, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, s0*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Activos Externos Netos (S^*_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax6);

sgtitle({'Choque Tecnológico Anticipado 5 Períodos (\chi = 0)',' '}, 'FontSize', 24, 'FontWeight', 'bold');
exportgraphics(fig2, 'FigB1.2.pdf')
fprintf('  -> Guardada: Figura_B1_Exp2_Tecnologico_Anticipado.png\n\n');

%% =========================================================================
% EXPERIMENTO 3: Choque tecnologico con tasa endogena (chi = 0.001) 
% =========================================================================
fprintf('Ejecutando Experimento 3: Choque tecnologico con tasa endogena ...\n');

set_param_value('chi', chi);
options_.periods = T;
M_.det_shocks = struct('exo_det', false, 'exo_id', 1, 'type', 'level', 'periods', 1:1, 'value', 5.0);
oo_ = perfect_foresight_setup(M_, options_, oo_);
[oo_, ~] = perfect_foresight_solver(M_, options_, oo_);

ct_3   = oo_.endo_simul(1, 2:Nplot+1)';
lt_3   = oo_.endo_simul(2, 2:Nplot+1)';
cat_3  = oo_.endo_simul(5, 2:Nplot+1)';
st_3   = oo_.endo_simul(6, 2:Nplot+1)';
rt_3   = oo_.endo_simul(7, 2:Nplot+1)';
tht_3  = oo_.endo_simul(8, 2:Nplot+1)';

fig3 = figure('Name', 'Tema B1.3', 'Position', fig_pos, 'Color', 'w');

ax1 = subplot(2, 3, 1);
plot(t_vec, tht_3, 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(t_vec, th0*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
plot(0, th0, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
xlim([0 50]);ylim("padded"); grid on;
title('Productividad (\theta_t)  [IMPULSO]', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.05 0.25 0.55]);
formatear_eje_impulso(ax1);
legend({'Trayectoria dinámica', 'Estado estacionario'}, 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');

ax2 = subplot(2, 3, 2);
plot(t_vec, ct_3, 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, css, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, css*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Consumo (C_t)', 'FontSize', 18, 'FontWeight', 'bold');
formatear_eje(ax2);

ax3 = subplot(2, 3, 3);
plot(t_vec, lt_3, 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, lss, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, lss*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Trabajo (L_t)', 'FontSize', 18, 'FontWeight', 'bold');
formatear_eje(ax3);

ax4 = subplot(2, 3, 4);
plot(t_vec, cat_3, 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, cass, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, cass*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Cuenta Corriente (CA_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax4);

ax5 = subplot(2, 3, 5);
plot(t_vec, st_3, 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, s0, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, s0*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Activos Externos Netos (S^*_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax5);

ax6 = subplot(2, 3, 6);
plot(t_vec, rt_3 * 100, 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, rlong*100, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, rlong*100*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Tasa de Interés Efectiva (r^*_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax6);

sgtitle({'Choque Tecnológico con Tasa de Interés Endógena (\chi = 0.001)',' '}, 'FontSize', 24, 'FontWeight', 'bold');
exportgraphics(fig3, 'FigB1.3.pdf')
fprintf('  -> Guardada: Figura_B1_Exp3_Tecnologico_TasaEndogena.png\n\n');

%% =========================================================================
% EXPERIMENTO 4: Choque a la tasa de interes no-anticipado 
% =========================================================================
fprintf('Ejecutando Experimento 4: Choque de tasa de interes no-anticipado ...\n');

rst_exp4 = rlong * ones(T, 1);
rst_exp4(1) = 0.10;
for t = 1:29
    rst_exp4(t+1) = rlong + 0.88 * (rst_exp4(t) - rlong);
end
tht_exp4 = th0 * ones(T, 1);

x0_exp4 = ones(T, 1) * [wss, lss, yss, cass, s0, css];
sol_exp4 = fsolve(@(x) ftrans_rbc(x, T, beta, phi, s0, rst_exp4, tht_exp4), x0_exp4, opts);

wt_4  = sol_exp4(:, 1);
lt_4  = sol_exp4(:, 2);
yt_4  = sol_exp4(:, 3);
cat_4 = sol_exp4(:, 4);
st_4  = sol_exp4(:, 5);
ct_4  = sol_exp4(:, 6);

fig4 = figure('Name', 'Tema B1.4', 'Position', fig_pos, 'Color', 'w');

ax1 = subplot(2, 3, 1);
plot(t_vec, rst_exp4(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(t_vec, rlong*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
plot(0, rlong, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
xlim([0 50]);ylim("padded"); grid on;
title('Tasa de Interés (r^*_t)  [IMPULSO]', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.05 0.25 0.55]);
formatear_eje_impulso(ax1);
legend({'Trayectoria dinámica', 'Estado estacionario'}, 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');

ax2 = subplot(2, 3, 2);
plot(t_vec, ct_4(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, css, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, css*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Consumo (C_t)', 'FontSize', 18, 'FontWeight', 'bold');
formatear_eje(ax2);

ax3 = subplot(2, 3, 3);
plot(t_vec, lt_4(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, lss, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, lss*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Trabajo (L_t)', 'FontSize', 18, 'FontWeight', 'bold');
formatear_eje(ax3);

ax4 = subplot(2, 3, 4);
plot(t_vec, yt_4(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, yss, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, yss*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Producto (Y_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax4);

ax5 = subplot(2, 3, 5);
plot(t_vec, cat_4(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, cass, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, cass*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Cuenta Corriente (CA_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax5);

ax6 = subplot(2, 3, 6);
plot(t_vec, st_4(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, s0, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, s0*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Activos Externos Netos (S^*_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax6);

sgtitle({'Choque a la Tasa de Interés No-Anticipado (\chi = 0)',' '}, 'FontSize', 24, 'FontWeight', 'bold');
exportgraphics(fig4, 'FigB1.4.pdf')
fprintf('  -> Guardada: Figura_B1_Exp4_TasaInteres_NoAnticipado.png\n\n');

%% =========================================================================
% EXPERIMENTO 5: Choque a la tasa de interes anticipado 
% =========================================================================
fprintf('Ejecutando Experimento 5: Choque de tasa de interes anticipado ...\n');

rst_exp5 = rlong * ones(T, 1);
rst_exp5(5) = 0.10;
for t = 5:34
    rst_exp5(t+1) = rlong + 0.88 * (rst_exp5(t) - rlong);
end
tht_exp5 = th0 * ones(T, 1);

x0_exp5 = ones(T, 1) * [wss, lss, yss, cass, s0, css];
sol_exp5 = fsolve(@(x) ftrans_rbc(x, T, beta, phi, s0, rst_exp5, tht_exp5), x0_exp5, opts);

wt_5  = sol_exp5(:, 1);
lt_5  = sol_exp5(:, 2);
yt_5  = sol_exp5(:, 3);
cat_5 = sol_exp5(:, 4);
st_5  = sol_exp5(:, 5);
ct_5  = sol_exp5(:, 6);

fig5 = figure('Name', 'Tema B1.5', 'Position', fig_pos, 'Color', 'w');

ax1 = subplot(2, 3, 1);
plot(t_vec, rst_exp5(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(t_vec, rlong*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
plot(0, rlong, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
xlim([0 50]);ylim("padded"); grid on;
title('Tasa de Interés (r^*_t)  [IMPULSO]', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.05 0.25 0.55]);
formatear_eje_impulso(ax1);
legend({'Trayectoria dinámica', 'Estado estacionario'}, 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');

ax2 = subplot(2, 3, 2);
plot(t_vec, ct_5(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, css, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, css*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Consumo (C_t)', 'FontSize', 18, 'FontWeight', 'bold');
formatear_eje(ax2);

ax3 = subplot(2, 3, 3);
plot(t_vec, lt_5(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, lss, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, lss*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Trabajo (L_t)', 'FontSize', 18, 'FontWeight', 'bold');
formatear_eje(ax3);

ax4 = subplot(2, 3, 4);
plot(t_vec, yt_5(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, yss, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, yss*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Producto (Y_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax4);

ax5 = subplot(2, 3, 5);
plot(t_vec, cat_5(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, cass, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, cass*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Cuenta Corriente (CA_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax5);

ax6 = subplot(2, 3, 6);
plot(t_vec, st_5(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, s0, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, s0*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Activos Externos Netos (S^*_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax6);

sgtitle({'Choque a la Tasa de Interés Anticipado 5 Períodos (\chi = 0)',' '}, 'FontSize', 24, 'FontWeight', 'bold');
exportgraphics(fig5, 'FigB1.5.pdf')
fprintf('  -> Guardada: Figura_B1_Exp5_TasaInteres_Anticipado.png\n\n');

%% =========================================================================
% EXPERIMENTO 6: Choque a la tasa de interes mundial con tasa endogena 
% =========================================================================
fprintf('Ejecutando Experimento 6: Choque de tasa con tasa endogena ...\n');

set_param_value('chi', chi);
options_.periods = T;
delta_rw = 0.10 - rlong;
M_.det_shocks = struct('exo_det', false, 'exo_id', 2, 'type', 'level', 'periods', 1:1, 'value', delta_rw);
oo_ = perfect_foresight_setup(M_, options_, oo_);
[oo_, ~] = perfect_foresight_solver(M_, options_, oo_);

ct_6   = oo_.endo_simul(1, 2:Nplot+1)';
lt_6   = oo_.endo_simul(2, 2:Nplot+1)';
cat_6  = oo_.endo_simul(5, 2:Nplot+1)';
st_6   = oo_.endo_simul(6, 2:Nplot+1)';
rt_6   = oo_.endo_simul(7, 2:Nplot+1)';
rwt_6  = oo_.endo_simul(9, 2:Nplot+1)';

fig6 = figure('Name', 'Tema B1.6', 'Position', fig_pos, 'Color', 'w');

ax1 = subplot(2, 3, 1);
plot(t_vec, rwt_6 * 100, 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(t_vec, rlong*100*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
plot(0, rlong*100, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
xlim([0 50]);ylim("padded"); grid on;
title('Tasa Mundial (r^w_t)  [IMPULSO]', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.05 0.25 0.55]);
formatear_eje_impulso(ax1);
legend({'Trayectoria dinámica', 'Estado estacionario'}, 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');

ax2 = subplot(2, 3, 2);
plot(t_vec, ct_6(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, css, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, css*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Consumo (C_t)', 'FontSize', 18, 'FontWeight', 'bold');
formatear_eje(ax2);

ax3 = subplot(2, 3, 3);
plot(t_vec, lt_6(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, lss, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, lss*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Trabajo (L_t)', 'FontSize', 18, 'FontWeight', 'bold');
formatear_eje(ax3);

ax4 = subplot(2, 3, 4);
plot(t_vec, cat_6(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, cass, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, cass*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Cuenta Corriente (CA_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax4);

ax5 = subplot(2, 3, 5);
plot(t_vec, st_6(1:Nplot), 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, s0, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, s0*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Activos Externos Netos (S^*_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax5);

ax6 = subplot(2, 3, 6);
plot(t_vec, rt_6 * 100, 'Color', col_dyn, 'LineWidth', lw_dyn); hold on;
plot(0, rlong*100, 'o', 'Color', col_ss, 'MarkerFaceColor', col_ss, 'MarkerSize', 5.5);
plot(t_vec, rlong*100*ones(Nplot, 1), '--', 'Color', col_ss, 'LineWidth', lw_ss);
xlim([0 50]);ylim("padded"); grid on;
title('Tasa de Interés Efectiva (r^*_t)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Períodos (t)', 'FontSize', 14, 'FontWeight', 'bold');
formatear_eje(ax6);

sgtitle({'Choque a la Tasa Mundial con Tasa Endógena (\chi = 0.001)',' '}, 'FontSize', 24, 'FontWeight', 'bold');
exportgraphics(fig6, 'FigB1.6.pdf')
fprintf('  -> Guardada: Figura_B1_Exp6_TasaInteres_TasaEndogena.png\n\n');

fprintf('=======================================================\n');
fprintf('  Todas las simulaciones se completaron con exito!\n');
fprintf('  Se generaron 6 figuras en formato 2x3.\n');
fprintf('=======================================================\n');

%% =========================================================================
% FUNCIONES AUXILIARES
% =========================================================================

function formatear_eje_impulso(ax)
    % Marco reforzado y destacado para el grafico del impulso
    set(ax, 'LineWidth', 1.0, 'Box', 'on', ...
            'FontSize', 14, 'TickDir', 'out', ...
            'XColor', [0.05 0.20 0.45], 'YColor', [0.05 0.20 0.45], ...
            'GridAlpha', 0.20, 'Color', [0.985 0.992 1.00]); % Fondo azul cielo
end

function formatear_eje(ax)
    % Marco estandar limpio para las variables de respuesta
    set(ax, 'LineWidth', 1.0, 'Box', 'on', ...
            'FontSize', 14, 'TickDir', 'out', ...
            'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2], ...
            'GridAlpha', 0.20, 'Color', 'w');
end

function F = ftrans_rbc(x, T, beta, phi, s0, rst, tht)
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
