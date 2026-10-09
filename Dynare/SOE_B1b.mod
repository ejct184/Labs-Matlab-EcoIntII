/*
 * SOE_B1b.mod
 *
 * Dynare counterpart of the small open economy WITH an endogenous, debt
 * elastic risk premium (chi = 0.001), which makes the model stationary.
 *
 * Experiments (numbered as in Matlab/TemaB1_experimentos.m):
 *   3 : unanticipated, persistent productivity shock (theta = 15 in t = 1)
 *   6 : unanticipated, persistent world interest rate shock (r* = 0.10 in t = 1)
 * Experiments 1, 2, 4 and 5 (chi = 0) are in SOE_B1.mod.
 *
 * TemaB1_experimentos.m solves experiments 3 and 6 with its own Dynare model
 * (temaB1.mod, not included in this repository). This file reproduces that
 * model exactly, written with the conventions of this repository:
 *   - Premium (eq. 8 of SOE.pdf): r(t) = rw(t) - chi*s(t), with s(t) the
 *     stock of foreign assets at the BEGINNING of period t (S(t-1) in
 *     temaB1.mod, where S is the end of period stock).
 *   - Shocks: AR(1) processes with persistence 0.88 hit in period 1 by an
 *     innovation (5 for productivity, 0.10 - r* for the world rate), with
 *     no truncation (unlike the fsolve experiments 1, 2, 4 and 5, where the
 *     shock lasts 30 periods). temaB1.mod solves the AR(1) equations as
 *     part of the model; here the same paths are built with the same
 *     recursion and passed as exogenous paths, which gives the same solution.
 *   - Terminal condition: Dynare's standard one (steady state in T+1), as in
 *     temaB1.mod. Running  dynare SOE_B1b -DT=1000  uses a longer horizon.
 * Note: the original Matlab/ftransB1b.m (experiment II of Matlab/TemaB1.m)
 * uses a different timing, r(t+1) = rst(t) - chi*s(t), the terminal
 * condition s(T-1) = s(T) and a truncated shock; see README.md.
 *
 * Solution method: Dynare perfect foresight (deterministic) solver over the
 * same T = 100 periods used in TemaB1_experimentos.m.
 *
 * Timing conventions:
 *   - Dynare period t is period t of the simulation; period 1 is the shock
 *     period (plotted as period 0).
 *   - ca(t) = y(t) - c(t) + r(t)*s(t) and s(t) = s(t-1) + ca(t-1); the end
 *     of period stock S(t) of temaB1.mod is s(t+1) = s(t) + ca(t), returned
 *     as st_end by get_TemaB1_paths.m.
 *   - Dynare period 0 (initial condition) is the steady state, where
 *     s = s0 = 0 and ca = 0, hence s(1) = s(0) + ca(0) = s0.
 *
 * Output: the struct soe_B1b with one field per experiment (soe_B1b.exp3,
 * soe_B1b.exp6), each with the paths wt, lt, yt, cat, st, st_end, ct, rste,
 * tht, rst (periods 1,...,T), the steady state ss_B1b, and the figures of
 * TemaB1_experimentos.m.
 *
 * Tested with Dynare 6.0 (Octave 8.4). Written for Dynare 5.x and 6.x.
 * Run with:  dynare SOE_B1b
 *            dynare SOE_B1b -DEXPERIMENTS=[6]     (only some experiments)
 */

@#ifndef T
@#define T = 100
@#endif
@#ifndef EXPERIMENTS
@#define EXPERIMENTS = [3, 6]
@#endif

var w  (long_name='Salario real')
    l  (long_name='Trabajo')
    y  (long_name='Producto')
    ca (long_name='Cuenta corriente')
    s  (long_name='Activos internacionales (inicio del periodo)')
    c  (long_name='Consumo')
    r  (long_name='Tasa de interes efectiva');

varexo theta (long_name='Productividad')
       rst   (long_name='Tasa de interes mundial');

parameters beta rlong phi s0 th0 chi;

// Calibration (identical to TemaB1_experimentos.m and temaB1.mod)
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

[name='fc: Euler equation']
c(+1)-beta*(1+r(+1))*c = 0;

[name='fr: effective interest rate with risk premium']
r = rst-chi*s;
end;

// Steady state of the model with risk premium (SOE.pdf, eqs. 15 and 20-22).
// With rst = rlong = 1/beta-1 it gives s = 0 = s0, i.e. exactly the initial
// steady state of the MATLAB codes and of temaB1.mod.
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
end;

steady;
resid;

// Initial conditions (period 0) are the steady state above. The MATLAB codes
// start from s(1) = s0, so the replication requires s0 to be the steady
// state level of assets (true in the MATLAB codes: s0 = 0 = s*).
if abs(oo_.steady_state(strcmp('s', M_.endo_names))-s0) > 1e-12, error('SOE_B1b: s0 must equal the steady state assets (rst-(1/beta-1))/chi, as in the MATLAB codes'); end

// Results of the experiments run in this call
soe_B1b = struct();

@#for IEXP in EXPERIMENTS
@#if IEXP != 3 && IEXP != 6
@#error "SOE_B1b.mod: only experiments 3 and 6 (experiments 1, 2, 4 and 5 are in SOE_B1.mod)"
@#endif

// ---------------------------------------------------------------------------
// Experiment @{IEXP}: AR(1) shock of temaB1.mod hit by one innovation in
// period 1 (TemaB1_experimentos.m sets M_.det_shocks for period 1 only)
// ---------------------------------------------------------------------------
theta_path = th0*ones(@{T},1);
rst_path   = rlong*ones(@{T},1);
@#if IEXP == 3
// Productivity: theta(1) = th0 + 5, then theta(t) - th0 = 0.88*(theta(t-1) - th0)
theta_path(1) = 15.0;
for t = 1:@{T-1}, theta_path(t+1) = th0+0.88*(theta_path(t)-th0); end
@#endif
@#if IEXP == 6
// World rate: rw(1) = rlong + (0.10 - rlong), then rw(t) - rlong = 0.88*(rw(t-1) - rlong)
rst_path(1) = 0.10;
for t = 1:@{T-1}, rst_path(t+1) = rlong+0.88*(rst_path(t)-rlong); end
@#endif

shocks(overwrite);
var theta;
periods 1:@{T};
values (theta_path);

var rst;
periods 1:@{T};
values (rst_path);
end;

perfect_foresight_setup(periods=@{T});
perfect_foresight_solver(tolf=1e-12, tolx=1e-12);

[soe_B1b.exp@{IEXP}, ss_B1b] = get_TemaB1_paths(M_, oo_);
plot_TemaB1(@{IEXP}, soe_B1b.exp@{IEXP}, ss_B1b);
@#endfor
