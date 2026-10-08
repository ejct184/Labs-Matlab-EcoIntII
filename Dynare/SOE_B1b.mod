/*
 * SOE_B1b.mod
 *
 * Dynare counterpart of Matlab/ftransB1b.m (experiment II of Matlab/TemaB1.m).
 *
 * Small open economy, infinite horizon, endogenous labor supply and an
 * ENDOGENOUS interest rate: a debt elastic risk premium (elasticity chi)
 * that makes the model stationary. The file computes the steady state and the
 * deterministic transition after a temporary, unanticipated increase in
 * productivity.
 *
 * Solution method: Dynare perfect foresight (deterministic) solver over the
 * same T = 100 periods used in TemaB1.m. With the default option
 * MATLAB_TERMINAL = 1, for every period t = 1,...,T the stacked Dynare system
 * is, equation by equation, the system F(x) = 0 coded in ftransB1b.m (700
 * equations in 700 unknowns).
 *
 * Timing conventions (identical to ftransB1b.m):
 *   - Dynare period t is row t of the fsolve solution in TemaB1.m.
 *     Period 1 is the shock period (it is plotted as period 0 in TemaB1.m).
 *   - s(t) is the stock of foreign assets at the BEGINNING of period t, so
 *     s(t) = s(t-1) + ca(t-1) and ca(t) = y(t) - c(t) + r(t)*s(t).
 *   - Interest rate: ftransB1b.m uses r(t+1) = rst(t) - chi*s(t) for t >= 1
 *     and r(1) = rst(1) (fr(1)). That is, r(t) = rst(t-1) - chi*s(t-1): the
 *     premium reacts to assets with a one period lag with respect to eq. (8)
 *     of SOE.pdf. The code (not the pdf) is replicated. The exogenous
 *     indicator d1 (= 1 only in period 1) imposes r(1) = rst(1) exactly.
 *   - Dynare period 0 (initial condition) is the steady state, where
 *     s = s0 = 0 and ca = 0, hence s(1) = s(0) + ca(0) = s0, as in
 *     fs(1) = st(1) - s0.
 *
 * Terminal condition: as ftransB1b.m, the last Euler equation is replaced by
 * s(T-1) = s(T) (fc(T) = st(T-1) - st(T)). It is imposed with the exogenous
 * indicator dT (= 1 only in period T). Because chi is small, the economy has
 * not fully returned to the steady state at T = 100, so this truncation
 * matters near T. Running
 *     dynare SOE_B1b -DMATLAB_TERMINAL=0 -DT=1000
 * instead uses Dynare's standard terminal condition (steady state at T+1)
 * with a long horizon, which approximates the infinite horizon solution.
 *
 * Tested with Dynare 6.0 (Octave 8.4). Written for Dynare 5.x and 6.x.
 * Run with:  dynare SOE_B1b
 */

@#ifndef T
@#define T = 100
@#endif
@#ifndef MATLAB_TERMINAL
@#define MATLAB_TERMINAL = 1
@#endif

var w  (long_name='Salario real')
    l  (long_name='Trabajo')
    y  (long_name='Producto')
    ca (long_name='Cuenta corriente')
    s  (long_name='Activos internacionales (inicio del periodo)')
    c  (long_name='Consumo')
    r  (long_name='Tasa de interes endogena');

varexo theta (long_name='Productividad')
       rst   (long_name='Tasa de interes mundial')
       d1    (long_name='Indicador del periodo inicial');
@#if MATLAB_TERMINAL
varexo dT    (long_name='Indicador del periodo terminal T');
@#endif

parameters beta rlong phi s0 th0 chi;

// Calibration (identical to TemaB1.m)
beta  = 0.95;           // discount factor
rlong = (1/beta)-1;     // long run real interest rate
phi   = 0.5;            // weight of leisure in the utility function
s0    = 0;              // initial foreign assets
th0   = 10;             // initial productivity
chi   = 0.001;          // elasticity of the risk premium

model;
[name='fw: wage equation']
w = theta;

[name='fl: labor supply']
l = 1-phi*(c/w);

[name='fy: production function']
y = theta*l;

[name='fca: current account']
ca = y-c+r*s;

[name='fs: law of motion of foreign assets']
s = s(-1)+ca(-1);

@#if MATLAB_TERMINAL
[name='fc: Euler equation (t<T) and terminal condition s(T-1)=s(T) (t=T)']
(1-dT)*(c(+1)-beta*(1+r(+1))*c) + dT*(s(-1)-s) = 0;
@#else
[name='fc: Euler equation']
c(+1)-beta*(1+r(+1))*c = 0;
@#endif

[name='fr: endogenous interest rate, r(1)=rst(1) (t=1)']
r = d1*rst + (1-d1)*(rst(-1)-chi*s(-1));
end;

// Steady state of the model with risk premium (SOE.pdf, eqs. 15 and 20-22).
// With rst = rlong = 1/beta-1 it gives s = 0 = s0, i.e. exactly the initial
// steady state used in TemaB1.m (css, wss, lss, yss, cass, s0, rlong).
steady_state_model;
r  = (1/beta)-1;                                     // Euler: 1 = beta*(1+r)
s  = (rst-r)/chi;                                    // r = rst - chi*s
c  = (th0/(1+phi))+(1/(1+phi))*r*s;                  // consumption
w  = th0;                                            // real wage
l  = 1-phi*(c/w);                                    // labor
y  = th0*l;                                          // output
ca = 0;                                              // current account
end;

initval;
theta = th0;
rst   = rlong;
d1    = 0;
@#if MATLAB_TERMINAL
dT    = 0;
@#endif
end;

steady;
resid;

// Initial conditions (period 0) are the steady state above. TemaB1.m starts
// experiment II from s(1) = s0, so the replication requires s0 to be the
// steady state level of assets (true in TemaB1.m: s0 = 0 = s*).
if abs(oo_.steady_state(strcmp('s', M_.endo_names))-s0) > 1e-12, error('SOE_B1b: s0 must equal the steady state assets (rst-(1/beta-1))/chi, as in TemaB1.m'); end

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

var d1;
periods 1;
values 1;

@#if MATLAB_TERMINAL
var dT;
periods @{T};
values 1;
@#endif
end;

perfect_foresight_setup(periods=@{T});
perfect_foresight_solver(tolf=1e-12, tolx=1e-12);

// Transition paths in the format of TemaB1.m (wt, lt, yt, cat, st, ct, rste,
// tht) and the same figure as figure(2) of TemaB1.m.
[soe_B1b, ss_B1b] = get_TemaB1_paths(M_, oo_);
plot_TemaB1(2, soe_B1b, ss_B1b);
