

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
Q = 1e4*eye(6);
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

%[Np_steps_H, Np_steps_h] = controllable_set(Hx, hx, Hu, hu, CIS_H, CIS_h, A_d, B_d,50);

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
    N=50

%% mpc e simulazione

T_sim=80

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
    

    % if isempty(delta_u_seq) || exitflag <= 0
    %     warning("Quadprog non ha trovato una soluzione ammissibile al passo %d (exitflag = %d)", tt, exitflag);
    %     delta_u_seq_first = zeros(3,1); % fallback
    % else
    %     delta_u_seq_first = delta_u_seq(1:3);
    % end

    if isempty(delta_u_seq) || exitflag <= 0
    warning("Quadprog fallito al passo %d (exitflag = %d), mantengo u precedente.", tt, exitflag);
    u_log(:, tt) = u_log(:, tt-1);  % ripeti l’ultimo controllo valido
    else
    delta_u_seq_first = delta_u_seq(1:3);
    u_log(:, tt) = u_ref + delta_u_seq_first;
    end

    %-------------PROBLEMA-----------------------------
    % delta_u_seq risulta vuoto, pertanto dice che vado ordine i limiti
    % dell'array...penso sia perchè quadprog non trova una soluzione
    % ottima.

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

    dxdt = @(t, x) dinamica_nonlineare(t, x, u_ref, k, k_ext, T_ext, C, tau);


    [~, xx] = ode45(dxdt, [0 Ts], x_current);

    x_log(:, tt+1) = xx(end, :)';

end

%% 5. Stampa risultati finali

x_final = x_log(:, end);  % Stato finale del sistema

% Separazione delle variabili
T1 = x_final(1);
T2 = x_final(2);
T3 = x_final(3);
Q1 = x_final(4);
Q2 = x_final(5);
Q3 = x_final(6);

fprintf('\n--- Risultati finali della simulazione ---\n');
fprintf('Temperature finali [°C]:\n');
fprintf('  T1 = %.2f\n', T1);
fprintf('  T2 = %.2f\n', T2);
fprintf('  T3 = %.2f\n', T3);

fprintf('\nPotenza termosifoni finali [W]:\n');
fprintf('  Q1 = %.2f\n', Q1);
fprintf('  Q2 = %.2f\n', Q2);
fprintf('  Q3 = %.2f\n', Q3);
fprintf('------------------------------------------\n');


%% 4. Plot dei risultati per il sistema termico

% Traslazione del CIS e del set N-step nelle coordinate originali
CIS_shifted = CIS_G + x_ref;  % x_ref è [289;289;289;100;100;100]
Np_step_set_shifted = Polyhedron(Np_steps_H, Np_steps_h) + x_ref;

% Proiezioni sui soli stati T1, T2, T3 (coordinate 1, 2, 3)
cis_temp = projection(CIS_shifted, [1 2 3]);
Np_step_temp = projection(Np_step_set_shifted, [1 2 3]);

figure
h_npstep = Np_step_temp.plot('Alpha', 0.05, 'LineWidth', 2, 'EdgeColor', 'blue');
hold on
h_cis = cis_temp.plot('Alpha', 0.1, 'EdgeColor', 'black');
h_traj = plot3(x_log(1, :), x_log(2, :), x_log(3, :), 'Color', [0 0 0.5]);
h_dots = scatter3(x_log(1,:), x_log(2,:), x_log(3,:), 30, 'cyan', 'filled');

title('Traiettoria del sistema termico')
xlabel('$T_1$ [$^\circ$C]', 'Interpreter', 'latex')
ylabel('$T_2$ [$^\circ$C]', 'Interpreter', 'latex')
zlabel('$T_3$ [$^\circ$C]', 'Interpreter', 'latex')

legend([h_cis, h_npstep, h_traj, h_dots], ...
    {'CIS', sprintf('%d-step set', N), 'Traiettoria', 'Campioni'}, ...
    'Interpreter','latex')

grid on
view(3)

% Plot anche per Q1 Q2 Q3
cis_power = projection(CIS_shifted, [4 5 6]);
Np_step_power = projection(Np_step_set_shifted, [4 5 6]);

figure
h_npstep = Np_step_power.plot('Alpha', 0.05, 'LineWidth', 2, 'EdgeColor', 'blue');
hold on
h_cis = cis_power.plot('Alpha', 0.1, 'EdgeColor', 'black');
h_traj = plot3(x_log(4, :), x_log(5, :), x_log(6, :), 'Color', [0.2 0 0.5]);
h_dots = scatter3(x_log(4,:), x_log(5,:), x_log(6,:), 30, 'magenta', 'filled');

title('Traiettoria delle potenze termosifoni')
xlabel('$Q_1$ [W]', 'Interpreter', 'latex')
ylabel('$Q_2$ [W]', 'Interpreter', 'latex')
zlabel('$Q_3$ [W]', 'Interpreter', 'latex')

legend([h_cis, h_npstep, h_traj, h_dots], ...
    {'CIS', sprintf('%d-step set', N), 'Traiettoria', 'Campioni'}, ...
    'Interpreter','latex')

grid on
view(3)




