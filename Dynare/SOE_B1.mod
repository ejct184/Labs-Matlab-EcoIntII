/*
 * SOE_B1.mod
 *
 * Dynare counterpart of Matlab/ftransB1.m (experiment I of Matlab/TemaB1.m).
 *
 * Small open economy, infinite horizon, endogenous labor supply and an
 * EXOGENOUS world interest rate (no endogenous risk premium). The file computes
 * the initial steady state and the deterministic transition after a temporary,
 * unanticipated increase in productivity.
 *
 * Solution method: Dynare perfect foresight (deterministic) solver over the
 * same T = 100 periods used in TemaB1.m. For every period t = 1,...,T the
 * stacked Dynare system is, equation by equation, the system F(x) = 0 coded
 * in ftransB1.m, so both programs solve exactly the same 600 equations in the
 * same 600 unknowns.
 *
 * Timing conventions (identical to ftransB1.m):
 *   - Dynare period t is row t of the fsolve solution in TemaB1.m.
 *     Period 1 is the shock period (it is plotted as period 0 in TemaB1.m).
 *   - s(t) is the stock of foreign assets at the BEGINNING of period t, so
 *     s(t) = s(t-1) + ca(t-1) and ca(t) = y(t) - c(t) + r*(t)*s(t).
 *   - Dynare period 0 (initial condition) is the pre-shock steady state,
 *     hence s(1) = s(0) + ca(0) = s0 + 0 = s0, as in fs(1) = st(1) - s0.
 *
 * Terminal condition: ftransB1.m replaces the last Euler equation by
 * s(T-1) = s(T) (fc(T) = st(T-1) - st(T)). This model has a unit root (any
 * level of assets is a steady state), so the new steady state depends on the
 * solution itself and cannot be fed to Dynare in advance through endval.
 * The same condition is therefore imposed with the exogenous indicator dT,
 * which equals 1 only in period T: in that period the Euler equation is
 * switched off and s(T-1) = s(T) is imposed instead. The terminal values of
 * c and rst (period T+1) are multiplied by zero and play no role.
 *
 * Tested with Dynare 6.0 (Octave 8.4). Written for Dynare 5.x and 6.x.
 * Run with:  dynare SOE_B1
 */

@#define T = 100

var w  (long_name='Salario real')
    l  (long_name='Trabajo')
    y  (long_name='Producto')
    ca (long_name='Cuenta corriente')
    s  (long_name='Activos internacionales (inicio del periodo)')
    c  (long_name='Consumo');

varexo theta (long_name='Productividad')
       rst   (long_name='Tasa de interes mundial')
       dT    (long_name='Indicador del periodo terminal T');

parameters beta rlong phi s0 th0;

// Calibration (identical to TemaB1.m)
beta  = 0.95;           // discount factor
rlong = (1/beta)-1;     // long run real interest rate
phi   = 0.5;            // weight of leisure in the utility function
s0    = 0;              // initial foreign assets
th0   = 10;             // initial productivity

model;
[name='fw: wage equation']
w = theta;

[name='fl: labor supply']
l = 1-phi*(c/w);

[name='fy: production function']
y = theta*l;

[name='fca: current account']
ca = y-c+rst*s;

[name='fs: law of motion of foreign assets']
s = s(-1)+ca(-1);

[name='fc: Euler equation (t<T) and terminal condition s(T-1)=s(T) (t=T)']
(1-dT)*(c(+1)-beta*(1+rst(+1))*c) + dT*(s(-1)-s) = 0;
end;

// Initial steady state (identical formulas to TemaB1.m)
steady_state_model;
c  = (th0/(1+phi))+(1/(1+phi))*((1-beta)/beta)*s0;   // consumption
w  = th0;                                            // real wage
l  = 1-phi*(c/w);                                    // labor
y  = th0*l;                                          // output
ca = 0;                                              // current account
s  = s0;                                             // foreign assets
end;

initval;
theta = th0;
rst   = rlong;
dT    = 0;
end;

steady;
resid;

// Temporary productivity increase, built exactly as in TemaB1.m:
// tht(1) = 15 and tht(i+1) = th0+0.88*(tht(i)-th0) for i = 1,...,29;
// productivity is back at th0 from period 31 onwards.
theta_path = th0*ones(@{T},1);
theta_path(1) = 15;
for i = 1:29, theta_path(i+1) = th0+0.88*(theta_path(i)-th0); end
theta_shock = theta_path(1:30);

shocks;
var theta;
periods 1:30;
values (theta_shock);

var dT;
periods @{T};
values 1;
end;

perfect_foresight_setup(periods=@{T});
perfect_foresight_solver(tolf=1e-12, tolx=1e-12);

// Transition paths in the format of TemaB1.m (wt, lt, yt, cat, st, ct, tht)
// and the same figure as figure(1) of TemaB1.m.
[soe_B1, ss_B1] = get_TemaB1_paths(M_, oo_);
plot_TemaB1(1, soe_B1, ss_B1);
