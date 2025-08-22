
clear; clc; close all

%% Richiamiamo lo script di inizializzazione
inizializzazione

%% Definizione delle matrici del costo quadratico
Np=30; %orizzonte di predizione
Tsim=20; %tempo di simulazione
Q = 1e3 * eye(6);
R = 1e1 * eye(3);
% S come soluzione di Riccati
[~, S] = dlqr(sys_discreto.A, sys_discreto.B, Q, R);

%% Verifica dell'esistenza del Controllable Invariant Set
[G, g] = cis2(sys_discreto.A, sys_discreto.B, zeros(6,1), zeros(3,1), Hx, hx, Hu, hu, Q, R);


%% Plot del CIS
CIS_G = Polyhedron(G, g);
CIS_G = minHRep(CIS_G);
CIS_G_T = projection(CIS_G, 1:3);
CIS_G_Q = projection(CIS_G, 4:6);

figure
subplot(1, 2, 1)
CIS_G_T.plot('Color', [0.2, 0.6, 0.8]); 
title("Proiezione del CIS delle temperature nelle stanze")
limitiTemp = [X_v_lin(7), X_v_lin(1)];
xlim(limitiTemp); ylim(limitiTemp); zlim(limitiTemp);
xlabel("T1 $[^{\circ}C]$", "Interpreter", "latex")
ylabel("T2 $[^{\circ}C]$", "Interpreter", "latex")
zlabel("T3 $[^{\circ}C]$", "Interpreter", "latex")

subplot(1, 2, 2)
CIS_G_Q.plot('Color', [0.2, 0.6, 0.8]); 
title("Proiezione del CIS della potenza termica dei termosifoni")
limitiQ = [X_v_lin(10), X_v_lin(4)];
xlim(limitiQ); ylim(limitiQ); zlim(limitiQ);
xlabel("Q1 $[W]$", "Interpreter", "latex")
ylabel("Q2 $[W]$", "Interpreter", "latex")
zlabel("Q3 $[W]$", "Interpreter", "latex")

%% N-step CS
x0_centrato = [284; 285; 284; 0; 10; 0] - x_ref(1:6);

%% Verifica fattibilità dal punto di partenza

%grafici sovrapposti
[Np_steps_H, Np_steps_h] = controllable_set2(Hx, hx, Hu, hu, G, g, A_d, B_d,Np);
Np_steps_set=Polyhedron('A',Np_steps_H,'b', Np_steps_h)
Np_step = Polyhedron(Np_steps_H , Np_steps_h);

%% Verifica fattibilità dal punto di partenza
Np_step = Polyhedron(Np_steps_H, Np_steps_h);
Np_step = Np_step.minHRep();
Np_steps_T = projection(Np_step , 1:3);
Np_steps_Q = projection(Np_step , 4:6);


figure
subplot(1 , 2 ,1)
Np_steps_T.plot('Color', [0.2, 0.6, 0.8]);
title("Temperature nelle stanze")
xlim(limitiTemp)
ylim(limitiTemp)
zlim(limitiTemp)
hold on
plot3(x0_centrato(1), x0_centrato(2), x0_centrato(3), '.', 'MarkerSize', 50, 'Color', 'r', 'DisplayName', 'Start Point')
xlabel("T1 $[^{\circ}C]$", "Interpreter", "latex")
ylabel("T2 $[^{\circ}C]$", "Interpreter", "latex")
zlabel("T3 $[^{\circ}C]$", "Interpreter", "latex")
legend_entries = {'N-Steps', 'Start Point'};
legend(legend_entries, 'Location', 'northeast', 'FontSize', 10)
hold off

subplot(1 , 2 ,2)
% Plot controllable set e salva handle
h_Np_Q = Np_steps_Q.plot('Color', [0.2, 0.6, 0.8]);
hold on

% Plot CIS e salva handle
h_CIS_Q = CIS_G_Q.plot();

title("Potenza termica dei termosifoni")
xlim(limitiQ)
ylim(limitiQ)
zlim(limitiQ)

% Trasparenze e colori usando gli handle
set(h_Np_Q, 'FaceColor', 'red', 'FaceAlpha', 0.3)
set(h_CIS_Q, 'FaceColor', 'blue', 'FaceAlpha', 0.3)

% Punto di partenza
plot3(x0_centrato(4), x0_centrato(5), x0_centrato(6), '.', 'MarkerSize', 50, 'Color', 'r', 'DisplayName', 'Start Point')

xlabel("Q1 $[W]$", "Interpreter", "latex")
ylabel("Q2 $[W]$", "Interpreter", "latex")
zlabel("Q3 $[W]$", "Interpreter", "latex")
legend_entries = {'N-Steps', 'Start Point'};
legend(legend_entries, 'Location', 'northeast', 'FontSize', 10)
hold off


%% Simulazione a tempo continuo con MPC

htt=[];
hxx = [];
u_online = [];
x_ini = [284 285 284 0 10 0]';
flag = zeros(1 , Np);
nsim=100;


[A_cal , A_cal_n , B_cal , B_cal_n,  Q_cal , R_cal] = Calligrafica(sys_discreto.A , sys_discreto.B , Q , R , S , Np);


for i = 1:nsim

    if i == 1
        x_run = x_ini-x_ref(1:6);
    else
        x_run = hxx(: , end)-x_ref(1:6);
    end

    [controlAction , flag(i)]= MPC(x_run , A_cal , A_cal_n , B_cal , B_cal_n,  Q_cal , R_cal , Np, G,g, X_v_lin, U_v_lin);
    tempo = linspace(Tsim*(i-1), Tsim*i , Tsim);
    controlAction = controlAction(1:3) + [100; 100; 100];
    u_online = [u_online,repmat(controlAction , 1 , Tsim)];
    dxdt = @(t,x) CasaEsterno(t, x, k, C, tau, T_ext, k_ext, controlAction);
    [tt, xx] = ode45(dxdt , tempo , x_run+x_ref);
    htt = [htt,tt'];
    hxx = [hxx,xx'];

end
%% Plot finale

tempo = htt/60; %[min]

% Configurazione globale per i plot
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultAxesFontWeight', 'bold');

% Prima figura: Evoluzioni degli stati
figure('Name', 'Evoluzioni degli stati', 'NumberTitle', 'off');
sgtitle("Evoluzioni degli stati nel tempo", 'FontSize', 16, 'FontWeight', 'bold');

% Subplot per le temperature
plot_T = subplot(2, 1, 1);
plot(tempo, hxx(1, :), 'r');  % rosso per T1
hold on;
plot(tempo, hxx(2, :), 'g');  % verde per T2
plot(tempo, hxx(3, :), 'b');  % blu per T3
yline(x_ref(1), 'Color', 'k', 'LineWidth', 0.5, 'Label', 'Obiettivo');  % nero per obiettivo, linea molto sottile
grid on;
legend(["T1 - Stanza 1", "T2 - Stanza 2", "T3 - Stanza 3", "Obiettivo"], ...
       'Location', 'best', 'FontSize', 11, 'FontWeight', 'bold');
ylabel("Temperatura $[^{\circ}C]$", "Interpreter", "latex", 'FontSize', 12);
xlabel("Tempo $[min]$", "Interpreter", "latex", 'FontSize', 12);
title("Evoluzione delle temperature nelle stanze", 'FontSize', 14, 'FontWeight', 'bold');
ylim([280, 290]);
xlim([0, max(tempo)]);
box on;

% Subplot per le potenze termiche
plot_Q = subplot(2, 1, 2);
plot(tempo, hxx(4, :), 'r');  % rosso per Q1
hold on;
plot(tempo, hxx(5, :), 'g');  % verde per Q2
plot(tempo, hxx(6, :), 'b');  % blu per Q3
yline(x_ref(4), 'Color', 'k', 'LineWidth', 0.5, 'Label', 'Obiettivo');  % nero per obiettivo, linea molto sottile
grid on;
legend(["Q1 - Termosifone 1", "Q2 - Termosifone 2", "Q3 - Termosifone 3", "Obiettivo"], ...
       'Location', 'best', 'FontSize', 11, 'FontWeight', 'bold');
ylabel("Potenza termica $[W]$", "Interpreter", "latex", 'FontSize', 12);
xlabel("Tempo $[min]$", "Interpreter", "latex", 'FontSize', 12);
title("Evoluzione della potenza termica dei termosifoni", 'FontSize', 14, 'FontWeight', 'bold');
ylim([0, 120]);
xlim([0, max(tempo)]);
box on;

% Seconda figura: Azioni di controllo
figure('Name', 'Azioni di controllo MPC', 'NumberTitle', 'off');
plot_U = plot(tempo, u_online(1, :), 'r');  % rosso per Q1
hold on;
plot(tempo, u_online(2, :), 'g');  % verde per Q2
plot(tempo, u_online(3, :), 'b');  % blu per Q3
grid on;
title("Azioni di controllo MPC nel tempo", 'FontSize', 16, 'FontWeight', 'bold');
ylabel("Potenza di controllo $[W]$", "Interpreter", "latex", 'FontSize', 12);
xlabel("Tempo $[min]$", "Interpreter", "latex", 'FontSize', 12);
legend(["Q1 - Controllo Stanza 1", "Q2 - Controllo Stanza 2", "Q3 - Controllo Stanza 3"], ...
       'Location', 'best', 'FontSize', 11, 'FontWeight', 'bold');
ylim([90, 150]);
xlim([0, max(tempo)]);
box on;

% Reset delle impostazioni globali
set(0, 'DefaultLineLineWidth', 'remove');
set(0, 'DefaultAxesFontSize', 'remove');
set(0, 'DefaultAxesFontWeight', 'remove');


