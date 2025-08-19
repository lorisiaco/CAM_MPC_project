

clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
set(0,'defaulttextInterpreter','latex')
rng('default');

%% Invariant set per il sistema controllato

addpath('funzioni');
main;

% Matrici del costo quadratico
Q = 1e3*eye(6);
R = 1e1*eye(3);                
% Control invariant set CIS_H*x <= CIS_h
[CIS_H, CIS_h] = cis(A_d, B_d, zeros(6,1), zeros(3,1), Hx, hx, Hu, hu, Q, R); 


CIS_G = Polyhedron(CIS_H, CIS_h);

% minimal H representation (elimina i vincoli ridondanti, ossia quelli che
% non  influenzano il set. Il risultato è equivalente ma con il minor
% numero di disuguaglianze lineari.
CIS_G = CIS_G.minHRep();
% crea una proiezione del poliedro CIS_G sul piano formato dalle variabili 
% 1 e 3 dello stato. Quindi, se lo stato x appartiene a R^6, sto
% proiettando il poliedro sui soli assi x1 e x3, ad esempio.
% RICORDA CHE: x = [T1 T2 T3 Q1 Q2 Q3]'
temperature = projection(CIS_G, [1,2,3]); %cubo T1 T2 T3
potenze = projection(CIS_G, [4,5,6]); % cubo Q1 Q2 Q3

dim = CIS_G.Dim;
disp(['dimensione poliedro: ', num2str(dim)]);

figure
subplot(1 , 2 , 1)
temperature.plot();
title("cis, temperatura stanze")
limT = [X_v_lin_min(1) X_v_lin_max(1)];
xlim(limT)
ylim(limT)
zlim(limT)
xlabel("T1 $[ ^{\circ}C]$" , Interpreter="latex")
ylabel("T2 $[ ^{\circ}C]$" , Interpreter="latex")
zlabel("T3 $[ ^{\circ}C]$" , Interpreter="latex")

subplot(1 , 2 , 2);
potenze.plot();
title("cis, potenza termosifoni")
limQ = [X_v_lin_min(5) X_v_lin_max(5)];
xlim(limQ)
ylim(limQ)
zlim(limQ)
xlabel("Q1 $[W]$" , Interpreter="latex")
ylabel("Q2 $[W]$" , Interpreter="latex")
zlabel("Q3 $[W]$" , Interpreter="latex")

% nota: ho cercato di traslare i grafici con le temperature e potenze reali
% (andando quindi ad aggiungere +273K e +100W; il problema però è che il
% modello realizzato è traslato...non penso si riesca a realizzare.

%% set controllabile in N passi


[Np_steps_H, Np_steps_h] = controllable_set(Hx, hx, Hu, hu, CIS_H, CIS_h, A_d, B_d,50);
Np_steps_set=Polyhedron('A',Np_steps_H,'b', Np_steps_h)
Np_step = Polyhedron(Np_steps_H , Np_steps_h);
Np_step = Np_step.minHRep();

Np_steps_T = projection(Np_step , 1:3);
Np_steps_Q = projection(Np_step , 4:6);


figure
subplot(1 , 2 ,1)
Np_steps_T.plot();
title("Temperature nelle stanze")
xlim(limT)
ylim(limT)
zlim(limT)
hold on
plot3(x0_centrato(1) ,x0_centrato(2), x0_centrato(3) , "." , MarkerSize=50)
xlabel("T1 $[ ^{\circ}C]$" , Interpreter="latex")
ylabel("T2 $[ ^{\circ}C]$" , Interpreter="latex")
zlabel("T3 $[ ^{\circ}C]$" , Interpreter="latex")

% legend(["n-steps" , "Punto di partenza"])

subplot(1 , 2 ,2)

% Plot controllable set e salva handle
h_Np_Q = Np_steps_Q.plot();
hold on

% Plot CIS e salva handle
h_CIS_Q = potenze.plot();

title("Potenza termica dei termosifoni")
xlim(limQ)
ylim(limQ)
zlim(limQ)

% Trasparenze e colori usando gli handle
set(h_Np_Q, 'FaceColor', 'red', 'FaceAlpha', 0.3)
set(h_CIS_Q, 'FaceColor', 'blue', 'FaceAlpha', 0.3)

% Punto di partenza
plot3(x0_centrato(1) ,x0_centrato(2), x0_centrato(3) , "." , MarkerSize=50)

xlabel("Q1 $[W]$" , Interpreter="latex")
ylabel("Q2 $[W]$" , Interpreter="latex")
zlabel("Q3 $[W]$" , Interpreter="latex")


%------------------------ PROBLEMA ----------------------------------

% problema: arrivati qui il programma si blocca senza terminare; 
% suppongo che il problema nasca dal fatto che fin dalla prima iterazione i
% vincoli sono veramente tanti (241!).
% Provando a cambiare parametri come N, Ts, R, Q il risultato rimane
% invariato...

%--------------------------------------------------------------------

% Intanto che cerco di sistemare questo problema mantengo 
    Np_steps_H=CIS_H
    Np_steps_h=CIS_h
    N=30

%% mpc e simulazione

T_sim=1

% Riferimento
x_ref = [289; 289; 289; 100; 100; 100];  % riferimento nel sistema reale
u_ref = [100; 100; 100];                 % potenza termosifoni a regime

% Riferimento del sistema linearizzato
x_ref_lin = zeros(6,1);  % origine, perché traslato x - x_ref
u_ref_lin = zeros(3,1);

mpc = mpc_ingredients(A_d,B_d,Hx,hx,Hu,hu,CIS_H,CIS_h,x_ref_lin,u_ref_lin,Q,R,N);

%matrice in cui salvo la traiettoria degli stati del sistema
% T_sim+1 perchè registro lo stato anche a t=0
x_log = zeros(6, T_sim+1);
%vettore dove salvo il controllo applicato ad ogni istant
u_log = zeros(3, T_sim);
% serve per salvare l'esito del quadProg
flags = zeros(1, T_sim);

%imposto lo stato iniziale del sistema
x_log(:, 1) = [284; 285; 284; 0; 10; 0];

for tt = 1:T_sim

    % prendo lo stato attuale linearizzato
    x_current = x_log(:, tt);

    x_lin = x_current - [289; 289; 289; 100; 100; 100];
    x_lin_shifted = x_lin - x_ref_lin;

    % Impostare i vincoli MPC in base alla condizione iniziale
    f = mpc.f_base * x_lin_shifted;
    b_ineq = mpc.b_ineq_base - mpc.b_ineq_x0_factor*x_lin_shifted;

    % risolve il problema di ottimizzazione quadratica
    [delta_u_seq, ~, exitflag] = quadprog(mpc.F, f, mpc.A_ineq, b_ineq);

    flags(tt) = exitflag;
    %trovata una sequenza ottimale di ingressi, prende solo il primo
    % PRINCIPIO RECEEDING HORIZON
    % in questo caso prendo i primi 3 poichè ho 3 ingressi
    delta_u_seq_first = delta_u_seq(1:3);

    %-------------PROBLEMA-----------------------------
    % delta_u_seq risulta vuoto, pertanto dice che vado ordine i limiti
    % dell'array...penso sia perchè quadprog non trova una soluzione
    % ottima. aumentando il numero di passi (20) però da una soluzione

    % il controllo è dato dal riferimento + la variazione ottima
    u_log(:,tt) = u_ref + delta_u_seq_first;


    %recap dei parametri passati in ingresso
    % u_current= valori di riferimento che il controllore da ai termosifoni.
    %(dai dati)
    % k= parametri di scambio termico tra stanze (dai dati) 
    % k_ext= parametro scambio termico tra stanze e esterno. (dai dati)
    % C= capacità termiche delle stanze (dai dati)
    % T_ext= temperatura esterna (dai dati)
    % tau= costanti di tempo dei termosifoni

    % dxdt = @(t, x) dinamica_nonlineare(t, x, u_ref, k, k_ext, T_ext, C, tau);
    % 
    % 
    % [~, xx] = ode45(dxdt, [0 Ts], x_current);
    % 
    % x_log(:, tt+1) = xx(end, :)';

end

