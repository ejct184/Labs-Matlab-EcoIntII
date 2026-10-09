// =========================================================================
// Tema B.1: Ciclos Economicos en una Economia Pequena y Abierta
// Curso: Economia Internacional II
// Profesor: Carlos Urrutia
// =========================================================================

// Variables endogenas:
// C     : Consumo
// L     : Oferta de trabajo
// Y     : Producto total
// w     : Salario real
// CA    : Cuenta corriente
// S     : Activos internacionales netos (stock a fin de periodo t)
// r     : Tasa de interes efectiva que enfrenta el pais (r*)
// theta : Choque de productividad (TFP)
// rw    : Tasa de interes real mundial exogena

var C L Y w CA S r theta rw;

// Variables exogenas (innovaciones a los procesos de choque):
varexo e_theta e_rw;

// Parametros:
parameters beta phi chi theta_bar r_bar rho_theta rho_rw;

beta      = 0.95;           // Factor de descuento intertemporal
phi       = 0.50;           // Peso del ocio en la funcion de utilidad
chi       = 0.001;          // Sensibilidad de la prima de riesgo al stock de deuda/activos
theta_bar = 10.0;           // Productividad en estado estacionario
r_bar     = (1/beta) - 1;   // Tasa de interes de estado estacionario (~5.263%)
rho_theta = 0.88;           // Persistencia del choque de productividad
rho_rw    = 0.88;           // Persistencia del choque a la tasa de interes

// Bloque de ecuaciones del modelo:
model;
    // 1. Salario real (determinacion competitiva: w_t = theta_t)
    w = theta;

    // 2. Oferta de trabajo (condicion estatica consumo-ocio)
    L = 1 - phi * (C / w);

    // 3. Funcion de produccion agregada
    Y = theta * L;

    // 4. Saldo de la cuenta corriente: CA_t = Y_t - C_t + r_t * S_{t-1}
    CA = Y - C + r * S(-1);

    // 5. Acumulacion de activos externos netos: S_t = S_{t-1} + CA_t
    // (En la convencion de Dynare, S_t se decide en t y es el stock al fin del periodo)
    S = S(-1) + CA;

    // 6. Ecuacion de Euler para el consumo intertemporal:
    // C_{t+1} = beta * (1 + r_{t+1}) * C_t
    C(+1) = beta * (1 + r(+1)) * C;

    // 7. Mecanismo de cierre: Tasa de interes endogena con prima por deuda
    // r_t = rw_t - chi * S_{t-1} (si chi=0, r_t = rw_t)
    r = rw - chi * S(-1);

    // 8. Proceso estocastico AR(1) para la productividad:
    theta - theta_bar = rho_theta * (theta(-1) - theta_bar) + e_theta;

    // 9. Proceso estocastico AR(1) para la tasa de interes mundial:
    rw - r_bar = rho_rw * (rw(-1) - r_bar) + e_rw;
end;

// Estado estacionario analitico:
initval;
    theta = theta_bar;
    rw    = r_bar;
    r     = r_bar;
    S     = 0;
    w     = theta_bar;
    L     = 1 / (1 + phi);
    Y     = theta_bar * L;
    C     = Y;
    CA    = 0;
    e_theta = 0;
    e_rw    = 0;
end;

steady;
check;
