
function F = ftransB1b(x)

global beta phi s0 chi
global T rst tht

% Secuencias a encontrar usando proceso iterativo

wt   = x(:,1);
lt   = x(:,2);
yt   = x(:,3);
cat  = x(:,4);
st   = x(:,5);
ct   = x(:,6);
rste = x(:,7);          % Tasa de interes endogena

% Inicializa sistema de ecuaciones

fw  = zeros(1,T);
fl  = zeros(1,T);
fy  = zeros(1,T);
fca = zeros(1,T);
fs  = zeros(1,T);  
fc  = zeros(1,T);  
fr  = zeros(1,T);  

% Ecuacion de salarios

for t=1:T
    fw(t) = wt(t)-tht(t);
end

% Oferta de trabajo

for t=1:T
    fl(t) = lt(t)-(1-phi*(ct(t)/wt(t)));
end

% Funcion de produccion

for t=1:T
    fy(t) = yt(t)-tht(t)*lt(t);
end

% Cuenta corriente

for t=1:T
    fca(t) = cat(t)-(yt(t)-ct(t)+rste(t)*st(t));
end

% Activos internacionales (dinamica)

for t=2:T 
    fs(t) = st(t)-(st(t-1)+cat(t-1));
end
fs(1)= st(1)-s0;

% Ecuacion de Euler (dinamica)

for t=1:T-1 
    fc(t) = ct(t+1)-beta*(1+rste(t+1))*ct(t);
end
fc(T)= st(T-1)-st(T);

% Ecuacion para la tasa de interes endogena

for t=1:T-1 
    fr(t+1) = rste(t+1)-(rst(t)-chi*st(t));
end
fr(1)= rste(1)-rst(1);

% Sistema final de ecuaciones

F = [fw;fl;fy;fca;fs;fc;fr];
F = F';

