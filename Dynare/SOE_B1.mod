/*
 * SOE_B1.mod
 *
 * Dynare counterpart of the small open economy WITHOUT risk premium
 * (chi = 0, exogenous world interest rate): Matlab/ftransB1.m, which is
 * identical to the function ftrans_rbc of Matlab/TemaB1_experimentos.m.
 *
 * Experiments (numbered as in Matlab/TemaB1_experimentos.m):
 *   1 : unanticipated, persistent productivity shock (theta = 15 in t = 1)
 *       (also experiment I of Matlab/TemaB1.m)
 *   2 : anticipated productivity shock (announced in t = 1, hits in t = 5)
 *   4 : unanticipated, persistent world interest rate shock (r* = 0.10 in t = 1)
 *   5 : anticipated world interest rate shock (announced in t = 1, hits in t = 5)
 * Experiments 3 and 6 (chi = 0.001) are in SOE_B1b.mod.
 *
 * Solution method: Dynare perfect foresight (deterministic) solver over the
 * same T = 100 periods used in the MATLAB codes. For every period
 * t = 1,...,T the stacked Dynare system is, equation by equation, the system
 * F(x) = 0 coded in ftransB1.m / ftrans_rbc, so both programs solve exactly
 * the same 600 equations in the same 600 unknowns. Each experiment only
 * changes the exogenous paths of productivity (theta) and of the world
 * interest rate (rst), which are built with the same code as in
 * TemaB1_experimentos.m and passed to Dynare for every period 1,...,T.
 * Anticipated shocks need nothing special: under perfect foresight the whole
 * path is known from period 1 on.
 *
 * Timing conventions (identical to the MATLAB codes):
 *   - Dynare period t is row t of the fsolve solution.
 *     Period 1 is the announcement/shock period (plotted as period 0).
 *   - s(t) is the stock of foreign assets at the BEGINNING of period t, so
 *     s(t) = s(t-1) + ca(t-1) and ca(t) = y(t) - c(t) + r*(t)*s(t).
 *   - Dynare period 0 (initial condition) is the pre-shock steady state,
 *     hence s(1) = s(0) + ca(0) = s0 + 0 = s0, as in fs(1) = st(1) - s0.
 *
 * Terminal condition: the MATLAB codes replace the last Euler equation by
 * s(T-1) = s(T) (fc(T) = st(T-1) - st(T)). This model has a unit root (any
 * level of assets is a steady state), so the new steady state depends on the
 * solution itself and cannot be fed to Dynare in advance through endval.
 * The same condition is therefore imposed with the exogenous indicator dT,
 * which equals 1 only in period T: in that period the Euler equation is
 * switched off and s(T-1) = s(T) is imposed instead. The terminal values of
 * c and rst (period T+1) are multiplied by zero and play no role.
 *
 * Output: the struct soe_B1 with one field per experiment (soe_B1.exp1,
 * soe_B1.exp2, soe_B1.exp4, soe_B1.exp5), each with the paths wt, lt, yt,
 * cat, st, ct, tht, rst (periods 1,...,T), the initial steady state ss_B1,
 * and the figures of TemaB1_experimentos.m.
 *
 * Tested with Dynare 6.0 (Octave 8.4). Written for Dynare 5.x and 6.x.
 * Run with:  dynare SOE_B1
 *            dynare SOE_B1 -DEXPERIMENTS=[2]     (only some experiments)
 */

@#ifndef T
@#define T = 100
@#endif
@#ifndef EXPERIMENTS
@#define EXPERIMENTS = [1, 2, 4, 5]
@#endif

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

// Calibration (identical to TemaB1.m and TemaB1_experimentos.m)
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

// Initial steady state (identical formulas to the MATLAB codes)
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

// Results of the experiments run in this call
soe_B1 = struct();

@#for IEXP in EXPERIMENTS
@#if IEXP != 1 && IEXP != 2 && IEXP != 4 && IEXP != 5
@#error "SOE_B1.mod: only experiments 1, 2, 4 and 5 (experiments 3 and 6 are in SOE_B1b.mod)"
@#endif

// ---------------------------------------------------------------------------
// Experiment @{IEXP}: exogenous paths built exactly as tht_exp@{IEXP} and
// rst_exp@{IEXP} in TemaB1_experimentos.m
// ---------------------------------------------------------------------------
theta_path = th0*ones(@{T},1);
rst_path   = rlong*ones(@{T},1);
@#if IEXP == 1
// Unanticipated persistent productivity shock
theta_path(1) = 15.0;
for t = 1:29, theta_path(t+1) = th0+0.88*(theta_path(t)-th0); end
@#endif
@#if IEXP == 2
// Anticipated productivity shock: known in t = 1, productivity rises in t = 5
theta_path(5) = 15.0;
for t = 5:34, theta_path(t+1) = th0+0.88*(theta_path(t)-th0); end
@#endif
@#if IEXP == 4
// Unanticipated persistent world interest rate shock
rst_path(1) = 0.10;
for t = 1:29, rst_path(t+1) = rlong+0.88*(rst_path(t)-rlong); end
@#endif
@#if IEXP == 5
// Anticipated world interest rate shock: known in t = 1, rate rises in t = 5
rst_path(5) = 0.10;
for t = 5:34, rst_path(t+1) = rlong+0.88*(rst_path(t)-rlong); end
@#endif

shocks(overwrite);
var theta;
periods 1:@{T};
values (theta_path);

var rst;
periods 1:@{T};
values (rst_path);

var dT;
periods @{T};
values 1;
end;

perfect_foresight_setup(periods=@{T});
perfect_foresight_solver(tolf=1e-12, tolx=1e-12);

[soe_B1.exp@{IEXP}, ss_B1] = get_TemaB1_paths(M_, oo_);
plot_TemaB1(@{IEXP}, soe_B1.exp@{IEXP}, ss_B1);
@#endfor
