clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
% set(0,'DefaultFigureWindowStyle', 'normal')
set(0,'defaulttextInterpreter','latex')
rng('default');

%% parametri modello

T_ext=273 %temperatura esterna costante [K]

% capacità termica stanza [J/s]
C1=6300 
C2=4600
C3=4200
C=[C1;C2;C3]


% parametro regolatore scambio di calore [W/K]
k1=16 
k2=18
k3=19
k=[k1;k2;k3]
k_ext=9

% la costante di tempo termosifone
tau1=580;     
tau2=520;
tau3=540;
tau=[tau1;tau2;tau3];

%stato iniziale

% la temperatura della i-esima stanza [K]
T1=284;     
T2=285;
T3=284;

% la potenza termica dell’i-esimo termosifone [W]
Q1=0; 
Q2=10;
Q3=0;

x_start=[T1 T2 T3 Q1 Q2 Q3]';


%stato obiettivo
x_ref=[289 289 289 100 100 100]'
u_ref=[100 100 100]'

%% linearizzo sistema 

syms T1 T2 T3 Q1 Q2 Q3 Q1r Q2r Q3r real

% real = ho a che fare con numeri reali
% comando syms = definisce simboli da usare nelle espressioni matematiche.

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

%rappresenta x_dot=f(x,u)
F = [f1; f2; f3; f4; f5; f6]; 
X = [T1 T2 T3 Q1 Q2 Q3];
U = [Q1r Q2r Q3r];

%derivata parziale rispetto agli stati
A_sym(X) = jacobian(F, X); 
%derivata parziale rispetto agli ingressi
B_sym(U) = jacobian(F, U);      

% Valutazione numerica delle matrici linearizzate in punto di equilibrio

A_lin = double(A_sym(x_ref(1), x_ref(2), x_ref(3), x_ref(4), x_ref(5), x_ref(6)));
B_lin = double(B_sym(u_ref(1), u_ref(2), u_ref(3)));

C_lin = eye(6);
D_lin = zeros(6, 3);

sys_lineare = ss(A_lin, B_lin, C_lin, D_lin);

x0_centrato = x_start - x_ref;

%% discretizzazione

% Parametri di discretizzazione
Ts = 1; 

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

%% vincoli

% Limiti fisici delle variabili (in unità assolute)
T_min = 282.5;   % Temperatura minima [K]
T_max = 300;     % Temperatura massima [K]

Q_min = 0;       % Potenza termosifoni minima (stati) [W]
Q_max = 150;     % Potenza termosifoni massima (stati) [W]

U_min = 0;       % Potenza di comando minima (ingressi) [W]
U_max = 150;     % Potenza di comando massima (ingressi) [W]

%% spiegazione

% Stiamo lavorando con variabili centrate nel punto d'equilibrio designato,
% per cui:
%    x_min - x_ref <= x_tilde <= x_max - x_ref
%    che equivale a: Hx * x_tilde <= hx
%
% Esempio: se T1 deve stare fra 282.5 e 300, ma il riferimento è 289,
% allora la temperatura centrata deve stare fra:
%    T1_tilde_min = 282.5 - 289 = -6.5
%    T1_tilde_max = 300   - 289 = +11
%
% Questa logica vale per tutte le variabili di stato:
%    - Le temperature (T1, T2, T3)
%    - Le potenze istantanee dei termosifoni (Q1, Q2, Q3)
%
% Idem per gli ingressi (Q1r, Q2r, Q3r): centrati rispetto a u_ref
%
% Per costruire i vincoli in forma matriciale:
%
%    Hx = [ eye(n);     --> per x_tilde <= x_max - x_ref
%           -eye(n) ];  --> per -x_tilde <= -(x_min - x_ref) => x_tilde >= x_min - x_ref
%
%    hx = [ x_max - x_ref;
%          -(x_min - x_ref) ];
%
% Lo stesso vale per Hu e hu per gli ingressi:
%    Hu * u_tilde <= hu
%
% Dove:
%    x_tilde = x - x_ref
%    u_tilde = u - u_ref
%
% In questo modo, anche se il sistema è linearizzato attorno a un punto non nullo,
% i vincoli rimangono corretti e coerenti nel dominio centrato usato dal controllo predittivo.


%% Vincoli sugli stati x = [T1; T2; T3; Q1; Q2; Q3]

% Costruiamo i vincoli assoluti sugli stati: prima i massimi, poi i minimi
X_max = [T_max; T_max; T_max; Q_max; Q_max; Q_max];  
X_min = [T_min; T_min; T_min; Q_min; Q_min; Q_min];  

% Centro i vincoli rispetto al punto di equilibrio (x_ref)
% Perché stiamo lavorando con deviazioni dal riferimento: x_tilde = x - x_ref

%forse l'errore risiede qui:

X_v_lin_max = X_max - x_ref;  % Vincoli superiori su x_tilde
X_v_lin_min = X_min - x_ref;  % Vincoli inferiori su x_tilde

% Metto insieme i vincoli in forma standard Hx * x <= hx
Hx = [ eye(6);     % x <= x_max
      -eye(6)];    % -x <= -x_min => x >= x_min

hx = [X_v_lin_max;
     -X_v_lin_min];   

%% Vincoli sugli ingressi u = [Q1r; Q2r; Q3r]

% Costruiamo i vincoli assoluti sugli ingressi
U_max_vec = [U_max; U_max; U_max];  
U_min_vec = [U_min; U_min; U_min];  

% Centro anche i vincoli sugli ingressi rispetto a u_ref
U_v_lin_max = U_max_vec - u_ref;  % Vincoli superiori su u_tilde
U_v_lin_min = U_min_vec - u_ref;  % Vincoli inferiori su u_tilde

% Vincoli in forma Hu * u <= hu
Hu = [ eye(3);     % u <= u_max
      -eye(3)];    % -u <= -u_min => u >= u_min

hu = [U_v_lin_max;
     -U_v_lin_min];  


