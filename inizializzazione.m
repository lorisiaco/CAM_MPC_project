clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
% set(0,'DefaultFigureWindowStyle', 'normal')
set(0,'defaulttextInterpreter','latex')
rng('default');
%% 0.PARAMETRI DELLA CASA

T_ext=273;  %la temperatura esterna, ipotizzata costante

C1=6300;     % la capacità termica della i-esima stanza
C2=4600;     
C3=4200;     
C=[C1;C2;C3];

k1=16;      % parametro che regola lo scambio di calore tra l’i-esima e la j-esima stanza
k2=18;
k3=19;
k=[k1;k2;k3];

k_ext=9;    % un parametro che regola lo scambio di calore tra le stanze della casa e l’esterno

tau1=580;     % la costante di tempo dell’i-esimo termosifone
tau2=520;
tau3=540;
tau=[tau1;tau2;tau3];

% Condizioni Iniziali

T1=284;     % la temperatura della i-esima stanza
T2=285;
T3=284;

Q1=0;       % la potenza termica dell’i-esimo termosifone
Q2=10;
Q3=0;

x_start=[T1 T2 T3 Q1 Q2 Q3]';

% Condizioni Obbiettivo 

x_ref=[289 289 289 100 100 100]';
u_ref=[100 100 100]';       %????? PERCHEEE ??????

%% 1.ODE SISTEMA
dxdt = @(t,x) CasaEsterno(t, x, k, C, tau, T_ext, k_ext, u_ref);

simulazione = 80 * 60;

% simulazione del comportamento del sistema
[tt, xx] = ode45(dxdt, linspace(0, simulazione, simulazione+1), x_start);
figure
hold on

tt = tt/60;

sgtitle("Simulazione del sistema con u di 100W") 

subplot(2,1,1)
plot(tt, xx(:,1), 'r', ...    % rosso per T1
     tt, xx(:,2), 'g', ...    % verde per T2
     tt, xx(:,3), 'b');       % blu per T3
title("Temperature delle stanze")
ylabel("Temperatura $[ ^{\circ}C]$" , Interpreter="latex")
xlabel("Tempo $[min]$",  Interpreter="latex")
legend(["T1" "T2" "T3"])

subplot(2,1,2)
plot(tt, xx(:,4), 'r', ...    % rosso per Q1
     tt, xx(:,5), 'g', ...    % verde per Q2
     tt, xx(:,6), 'b');       % blu per Q3
title("Potenza dei termosifoni")
ylabel("Potenza $[W]$" , Interpreter="latex")
xlabel("Tempo $[min]$",  Interpreter="latex")
legend(["Q1" "Q2" "Q3"])

hold off

%% 2.LINEARIZZAZIONE SISTEMA

syms T1 T2 T3 Q1 Q2 Q3 Q1r Q2r Q3r real

% Scambi termici non lineari tra stanze
k12 = k(1) + 4 / (1 + exp(-0.5 * abs(T1 - T2)));
k13 = k(1) + 4 / (1 + exp(-0.5 * abs(T1 - T3)));
k23 = k(2) + 4 / (1 + exp(-0.5 * abs(T2 - T3)));

% Equazioni dinamiche
f1 = (Q1 - k12 * (T1 - T2) - k13 * (T1 - T3) - k_ext * (T1 - T_ext)) / C(1);
f2 = (Q2 + k12 * (T1 - T2) - k23 * (T2 - T3) - k_ext * (T2 - T_ext)) / C(2);
f3 = (Q3 + k13 * (T1 - T3) + k23 * (T2 - T3) - k_ext * (T3 - T_ext)) / C(3);

% Dinamica dei termosifoni
f4 = (Q1r - Q1) / tau(1);
f5 = (Q2r - Q2) / tau(2);
f6 = (Q3r - Q3) / tau(3);

F = [f1; f2; f3; f4; f5; f6];
X = [T1 T2 T3 Q1 Q2 Q3];
U = [Q1r Q2r Q3r];

A_sym(X) = jacobian(F, X);      % derivata di F rispetto agli stati
B_sym(U) = jacobian(F, U);      % derivata di F rispetto agli ingressi

% Valutazione numerica delle matrici linearizzate in punto di equilibrio
A_lin = double(A_sym(x_ref(1), x_ref(2), x_ref(3), x_ref(4), x_ref(5), x_ref(6)));
B_lin = double(B_sym(u_ref(1), u_ref(2), u_ref(3)));

C_lin = eye(6);
D_lin = zeros(6, 3);

sys_lineare = ss(A_lin, B_lin, C_lin, D_lin);

x0_centrato = x_start - x_ref;

disp("Autovalori della matrice A linearizzata:")
disp(eig(A_lin))


%% 3.DISCRETIZZAZIONE DEL SISTEMA

% Parametri di discretizzazione
Ts = 60;  % Periodo di campionamento in secondi (1 minuto)

% Discretizzazione del sistema lineare usando metodo zero-order hold
sys_discreto = c2d(sys_lineare, Ts, 'zoh');

% Estrazione delle matrici del sistema discreto
A_d = sys_discreto.A;
B_d = sys_discreto.B;
C_d = sys_discreto.C;
D_d = sys_discreto.D;

% Verifica della stabilità del sistema discreto
autovalori_discreti = eig(A_d);
disp("Autovalori della matrice A discreta:")
disp(autovalori_discreti)

%% 4.MATRICE RAGGIUNGIBILITA' 

% Analisi della matrice di controllabilità per il sistema lineare continuo
n = size(A_lin, 1);  % Dimensione dello stato
m = size(B_lin, 2);  % Dimensione dell'ingresso

% Matrice di controllabilità per il sistema continuo 
Mr_c = ctrb(A_lin, B_lin);

% Rango della matrice di controllabilità continua
R_c = rank(Mr_c);
disp("Analisi Controllabilità Sistema Continuo:")
disp("Dimensione stato: " + num2str(n))
disp("Dimensione ingresso: " + num2str(m))
disp("Rango matrice controllabilità: " + num2str(R_c))
if R_c == n
    disp("Il sistema continuo è completamente controllabile")
else
    disp("Il sistema continuo NON è completamente controllabile")
    disp("   Stati non controllabili: " + num2str(n - R_c))
end
disp(' ');

% Analisi della matrice di controllabilità per il sistema discreto
% Matrice di controllabilità per il sistema discreto usando ctrb
Mr_d = ctrb(A_d, B_d);

% Rango della matrice di controllabilità discreta
R_d = rank(Mr_d);
disp(' ');
disp("Analisi Controllabilità Sistema Discreto:")
disp("Rango matrice controllabilità: " + num2str(R_d))
if R_d == n
    disp("Il sistema discreto è completamente controllabile")
else
    disp("Il sistema discreto NON è completamente controllabile")
    disp("   Stati non controllabili: " + num2str(n - R_d))
end

% Confronto tra le due matrici di controllabilità
disp(' ');
disp("Confronto Matrici di Controllabilità:")
disp("Rango continuo: " + num2str(R_c))
disp("Rango discreto: " + num2str(R_d))

if R_c == R_d
    disp("La discretizzazione mantiene la controllabilità")
else
    disp("ATTENZIONE: La discretizzazione ha modificato la controllabilità")
end

%% 5.VINCOLI

% Vincoli fisici del sistema
T_v = [300, 282.5];    % [T_max, T_min] - Limiti temperatura (°K)
Q_v = [150, 0];        % [Q_max, Q_min] - Limiti potenza termosifoni (W)
U_v = [150, 0];        % [U_max, U_min] - Limiti potenza di riferimento (W)

% Creazione vincoli per stati e ingressi
% Vincoli stati: [T_max, T_max, T_max, Q_max, Q_max, Q_max, T_min, T_min, T_min, Q_min, Q_min, Q_min]
X_v = [T_v(1) * ones(3,1);    % Limiti superiori temperature
          Q_v(1) * ones(3,1);     % Limiti superiori potenze
          T_v(2) * ones(3,1);     % Limiti inferiori temperature  
          Q_v(2) * ones(3,1)];    % Limiti inferiori potenze

% Vincoli ingressi: [U_max, U_max, U_max, U_min, U_min, U_min]
U_v = [U_v(1) * ones(3,1);    % Limiti superiori ingressi
          U_v(2) * ones(3,1)];    % Limiti inferiori ingressi

% Centratura dei vincoli nel punto di equilibrio
X_v_lin = X_v - [x_ref; x_ref];  % Vincoli stati centrati
U_v_lin = U_v - [u_ref; u_ref];  % Vincoli ingressi centrati

% Creazione matrici vincolo per MPC
Hx = [eye(6); -eye(6)];  % Matrice vincoli stati
hx = [ones(6,1); -ones(6,1)] .* X_v_lin;  % Vettore vincoli stati

Hu = [eye(3); -eye(3)];  % Matrice vincoli ingressi  
hu = [ones(3,1); -ones(3,1)] .* U_v_lin;  % Vettore vincoli ingressi





