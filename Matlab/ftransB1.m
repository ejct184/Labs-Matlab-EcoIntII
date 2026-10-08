
function F = ftransB1(x)

global beta phi s0
global T rst tht

% Secuencias a encontrar usando proceso iterativo

wt  = x(:,1);
lt  = x(:,2);
yt  = x(:,3);
cat = x(:,4);
st  = x(:,5);
ct  = x(:,6);

% Inicializa sistema de ecuaciones

fw  = zeros(1,T);
fl  = zeros(1,T);
fy  = zeros(1,T);
fca = zeros(1,T);
fs  = zeros(1,T);  
fc  = zeros(1,T);  

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
    fca(t) = cat(t)-(yt(t)-ct(t)+rst(t)*st(t));
end

% Activos internacionales (dinamica)

for t=2:T 
    fs(t) = st(t)-(st(t-1)+cat(t-1));
end
fs(1)= st(1)-s0;

% Ecuacion de Euler (dinamica)

for t=1:T-1 
    fc(t) = ct(t+1)-beta*(1+rst(t+1))*ct(t);
end
fc(T)= st(T-1)-st(T);

% Sistema final de ecuaciones

F = [fw;fl;fy;fca;fs;fc]; %6x100
F = F'; %100x6

