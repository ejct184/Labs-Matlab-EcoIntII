
% Resuelve el estado estacionario y la transicion del modelo de horizonte
% infinito con oferta de trabajo endogena frente a choques no anticipados
% de productividad (el codigo puede adaptarse para los choques anticipados
% y para los cambios en la tasa de interes).
%
% Tambien produce algunos graficos del Tema B1 visto en clase.
%
% El codigo llama a la funcion ftransB1 con el sistema de ecuaciones que
% caracteriza la dinamica del modelo. La funcion ftransB1b añade a este
% sistema de ecuaciones el ajuste de la tasa de interes mundial endogena
% que asegura la estacionareidad del modelo.

clear

global beta phi s0 chi
global T rst tht

% Lista de parametros

beta  = 0.95;           % Factor de descuento
rlong = (1/beta)-1;     % Tasa de inetrés real de largo plazo
phi   = 0.5;            % Peso del ocio en la funcion de utilidad

s0    = 0;              % Activos internacionales iniciales
th0   = 10;             % Productividad inicial

options = optimset('Display','off');
                              
% Resuelve estado estacionario inicial

css   = (th0/(1+phi))+(1/(1+phi))*((1-beta)/beta)*s0;   % consumo
wss   = th0;                                            % salario real
lss   = 1-phi*(css/wss);                                % trabajo
yss   = th0*lss;                                        % producto
cass  = 0;                                          % cuenta corriente

%% I) Resuelve experimento: Aumento temporal en productividad

T  = 100;                     % Periodos para llegar al estado estacionario

rst = rlong*ones(T,1);        % Inicializa los dos choques (tasa de interes
tht = th0*ones(T,1);          % y productividad)

tht(1) = 15;                  % Introduce aumento temporal en productividad
for i = 1:29
    tht(i+1)=th0+0.88*(tht(i)-th0);
end

x0 = ones(1,T)'*[wss;lss;yss;cass;s0;css]'; %100x1 1x6

[y,Fval,exitflag] = fsolve(@ftransB1,x0,options);

disp(' ')
if exitflag==1
    disp('Algoritmo convergio a la solucion')
else
    disp('Algoritmo no convergio correctamente')    
end

wt  = y(:,1);
lt  = y(:,2);
yt  = y(:,3);
cat = y(:,4);
st  = y(:,5);
ct  = y(:,6);

% Grafica resultados

figure(1)

subplot(3,2,1)
plot((1:51)-1,tht(1:51),'b-',0,th0,'ro',(1:51)-1,th0*ones(51,1),'r--')
axis ([0 50 8 17])
xlabel('Productividad')

subplot(3,2,2)
plot((1:51)-1,ct(1:51),'b-',0,css,'ro',(1:51)-1,css*ones(51,1),'r--')
axis ([0 50 5 9])
xlabel('Consumo')

subplot(3,2,3)
plot((1:51)-1,lt(1:51),'b-',0,lss,'ro',(1:51)-1,lss*ones(51,1),'r--')
axis ([0 50 0.6 0.8])
xlabel('Trabajo')
subplot(3,2,4)

plot((1:51)-1,yt(1:51),'b-',0,yss,'ro',(1:51)-1,yss*ones(51,1),'r--')
axis ([0 50 5 12])
xlabel('Producto')

subplot(3,2,5)
plot((1:51)-1,cat(1:51),'b-',0,cass,'ro',(1:51)-1,cass*ones(51,1),'r--')
axis ([0 50 -2 5])
xlabel('Cuenta Corriente')

subplot(3,2,6)
plot((1:51)-1,st(1:51),'b-',0,s0,'ro',(1:51)-1,s0*ones(51,1),'r--')
axis ([0 50 -10 35])
xlabel('Activos')

%% II) Repite experimento con tasa de interes endogena

chi = 0.001;        % Elasticidad de la prima de riesgo

T  = 100;

rst = rlong*ones(T,1);
tht = th0*ones(T,1);

tht(1) = 15;
for i = 1:29
    tht(i+1)=th0+0.88*(tht(i)-th0);
end

x0 = ones(1,T)'*[wss;lss;yss;cass;s0;css;rlong]';

[y,Fval,exitflag] = fsolve(@ftransB1b,x0,options);

disp(' ')
if exitflag==1
    disp('Algoritmo convergio a la solucion')
else
    disp('Algoritmo no convergio correctamente')    
end
disp(' ')

wt  = y(:,1);
lt  = y(:,2);
yt  = y(:,3);
cat = y(:,4);
st  = y(:,5);
ct  = y(:,6);
rste = y(:,7);

% Grafica resultados

figure(2)

subplot(3,2,1)
plot((1:51)-1,tht(1:51),'b-',0,th0,'ro',(1:51)-1,th0*ones(51,1),'r--')
axis ([0 50 8 17])
xlabel('Productividad')

subplot(3,2,2)
plot((1:51)-1,ct(1:51),'b-',0,css,'ro',(1:51)-1,css*ones(51,1),'r--')
axis ([0 50 5 9])
xlabel('Consumo')

subplot(3,2,3)
plot((1:51)-1,lt(1:51),'b-',0,lss,'ro',(1:51)-1,lss*ones(51,1),'r--')
axis ([0 50 0.6 0.8])
xlabel('Trabajo')

subplot(3,2,4)
plot((1:51)-1,cat(1:51),'b-',0,cass,'ro',(1:51)-1,cass*ones(51,1),'r--')
axis ([0 50 -2 5])
xlabel('Cuenta Corriente')

subplot(3,2,5)
plot((1:51)-1,st(1:51),'b-',0,s0,'ro',(1:51)-1,s0*ones(51,1),'r--')
axis ([0 50 -10 35])
xlabel('Activos')

subplot(3,2,6)
plot((1:51)-1,rste(1:51)*100,'b-',0,rlong*100,'ro',(1:51)-1,rlong*100*ones(51,1),'r--')
axis ([0 50 3 12])
xlabel('Tasa de Interes')
