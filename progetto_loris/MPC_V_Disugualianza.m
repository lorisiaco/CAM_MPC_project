% In questo script valutiamo MPC con vincolo di disuguaglianza, quindi con un control invariant set (CIS)
% * Se Q > 1e2 probabilmente il calcolo del N Steps controllable set si
% blocca quindi in quel caso saltare la parte e impostare manualmente i
% passi facendo: Np = #passi

clear; clc; close all

%% Impostazioni dello script
Ts = 60; % [secondi] - tempo di campionamento

%% Richiamiamo lo script di inizializzazione
inizializzazione

%% Definizione delle matrici del costo quadratico
Q = 1.e3 * eye(6);
R = 1e1 * eye(3);
% S come soluzione di Riccati
[~, S] = dlqr(sys_discreto.A, sys_discreto.B, Q, R);

%% Verifica dell'esistenza del Controllable Invariant Set
[G, g] = CIS(sys_discreto.A, sys_discreto.B, zeros(6,1), zeros(3,1), Hx, hx, Hu, hu, Q, R);

%% Plot del CIS
CIS_G = Polyhedron(G, g);
CIS_G = minHRep(CIS_G);
CIS_G_T = projection(CIS_G, 1:3);
CIS_G_Q = projection(CIS_G, 4:6);

figure
subplot(1, 2, 1)
CIS_G_T.plot('Color', [0.2, 0.6, 0.8]); % Colore blu elegante invece del rosso
title("Proiezione del CIS delle temperature nelle stanze")
limitiTemp = [X_v_lin(7), X_v_lin(1)];
xlim(limitiTemp); ylim(limitiTemp); zlim(limitiTemp);
xlabel("T1 $[^{\circ}C]$", "Interpreter", "latex")
ylabel("T2 $[^{\circ}C]$", "Interpreter", "latex")
zlabel("T3 $[^{\circ}C]$", "Interpreter", "latex")

subplot(1, 2, 2)
CIS_G_Q.plot('Color', [0.2, 0.6, 0.8]); % Colore blu elegante invece del rosso
title("Proiezione del CIS della potenza termica dei termosifoni")
limitiQ = [X_v_lin(10), X_v_lin(4)];
xlim(limitiQ); ylim(limitiQ); zlim(limitiQ);
xlabel("Q1 $[W]$", "Interpreter", "latex")
ylabel("Q2 $[W]$", "Interpreter", "latex")
zlabel("Q3 $[W]$", "Interpreter", "latex")

%% N-step controllable set
x0_centrato = [284; 285; 284; 0; 10; 0] - x_ref(1:6);

% [Np_steps_H, Np_steps_h, Np] = CS(Hx, hx, Hu, hu, G, g, sys_discreto.A, sys_discreto.B, x0_centrato);
Np = 10;
Np_steps_H = G;
Np_steps_h = g;

disp(" ")
disp("Passi minimi per entrare nel CIS: " + Np);

%% Verifica fattibilità dal punto di partenza
trasp = 0.3; %togliere serve per le fugure
Np_step = Polyhedron(Np_steps_H, Np_steps_h);
Np_step = Np_step.minHRep();
Np_steps_T = projection(Np_step, 1:3);
Np_steps_Q = projection(Np_step, 4:6);

figure
subplot(1, 2, 1)
h1 = Np_steps_T.plot('Color', [0.2, 0.6, 0.8]); % Colore blu elegante
alpha(h1, 0.3); % Imposta trasparenza
title("Temperature nelle stanze")
xlim(limitiTemp); ylim(limitiTemp); zlim(limitiTemp);
hold on
% Punto di partenza
plot3(x0_centrato(1), x0_centrato(2), x0_centrato(3), '.', 'MarkerSize', 50, 'Color', 'r', 'DisplayName', 'Start Point')
% Traiettoria simulata (se disponibile)
if exist('hxx', 'var') && ~isempty(hxx)
    % Plot della traiettoria delle temperature
    traj_temp = hxx(1:3, :) - x_ref(1:3);
    plot3(traj_temp(1, :), traj_temp(2, :), traj_temp(3, :), 'k-', 'LineWidth', 2, 'DisplayName', 'Traiettoria simulata')
    % Punto finale
    plot3(traj_temp(1, end), traj_temp(2, end), traj_temp(3, end), 'go', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'DisplayName', 'Punto finale')
end
xlabel("T1 $[^{\circ}C]$", "Interpreter", "latex")
ylabel("T2 $[^{\circ}C]$", "Interpreter", "latex")
zlabel("T3 $[^{\circ}C]$", "Interpreter", "latex")
% Legenda personalizzata per evitare "data1"
legend_entries = {'N-Steps', 'Start Point'};
if exist('hxx', 'var') && ~isempty(hxx)
    legend_entries = [legend_entries, {'Traiettoria simulata', 'Punto finale'}];
end
legend(legend_entries, 'Location', 'northeast', 'FontSize', 10)
hold off

subplot(1, 2, 2)
h2 = Np_steps_Q.plot('Color', [0.2, 0.6, 0.8]); % Colore blu elegante
alpha(h2, 0.3); % Imposta trasparenza
title("Potenza termica dei termosifoni")
xlim(limitiQ); ylim(limitiQ); zlim(limitiQ);
hold on
% Punto di partenza
plot3(x0_centrato(4), x0_centrato(5), x0_centrato(6), '.', 'MarkerSize', 50, 'Color', 'r', 'DisplayName', 'Start Point')
% Traiettoria simulata (se disponibile)
if exist('hxx', 'var') && ~isempty(hxx)
    % Plot della traiettoria delle potenze
    traj_pot = hxx(4:6, :) - x_ref(4:6);
    plot3(traj_pot(1, :), traj_pot(2, :), traj_pot(3, :), 'k-', 'LineWidth', 2, 'DisplayName', 'Traiettoria simulata')
    % Punto finale
    plot3(traj_pot(1, end), traj_pot(2, end), traj_pot(3, end), 'go', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'DisplayName', 'Punto finale')
end
xlabel("Q1 $[W]$", "Interpreter", "latex")
ylabel("Q2 $[W]$", "Interpreter", "latex")
zlabel("Q3 $[W]$", "Interpreter", "latex")
% Legenda personalizzata per evitare "data1"
legend_entries = {'N-Steps', 'Start Point'};
if exist('hxx', 'var') && ~isempty(hxx)
    legend_entries = [legend_entries, {'Traiettoria simulata', 'Punto finale'}];
end
legend(legend_entries, 'Location', 'northeast', 'FontSize', 10)
hold off

%% Simulazione a tempo continuo con MPC
if ~exist("Np" , "var")
    Np = 10;
end

htt=[];
hxx = [];
u_online = [];
x_ini = [284 285 284 0 10 0]';
flag = zeros(1 , Np);
n_sim = 100;

[A_cal , A_cal_n , B_cal , B_cal_n,  Q_cal , R_cal] = Calligrafica(sys_discreto.A , sys_discreto.B , Q , R , S , Np);


for i = 1:n_sim

    if i == 1
        x_run = x_ini-x_ref(1:6);
    else
        x_run = hxx(: , end)-x_ref(1:6);
    end

    [controlAction , flag(i)]= MPC(x_run , A_cal , A_cal_n , B_cal , B_cal_n,  Q_cal , R_cal , Np, G,g, X_v_lin, U_v_lin);
    tempo = linspace(Ts*(i-1), Ts*i , Ts);
    controlAction = controlAction(1:3) + [100; 100; 100];
    u_online = [u_online,repmat(controlAction , 1 , Ts)];
    dxdt = @(t,x) CasaEsterno(t, x, k, C, tau, T_ext, k_ext, controlAction);
    [tt, xx] = ode45(dxdt , tempo , x_run+x_ref);
    htt = [htt,tt'];
    hxx = [hxx,xx'];

end
%% Plot finale

% Plot diretto senza utilizzare la funzione plotSimulazione
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


%------------------------------------------------------------------
%grafico in più


% % Traslazione del CIS e del set N-step nelle coordinate originali
% CIS_shifted = CIS_G + x_ref;  % x_ref è [289;289;289;100;100;100]
% Np_step_set_shifted = Polyhedron(Np_steps_H, Np_steps_h) + x_ref;
% 
% % Proiezioni sui soli stati T1, T2, T3 (coordinate 1, 2, 3)
% cis_temp = projection(CIS_shifted, [1 2 3]);
% Np_step_temp = projection(Np_step_set_shifted, [1 2 3]);
% 
% figure
% h_npstep = Np_step_temp.plot('Alpha', 0.05, 'LineWidth', 2, 'EdgeColor', 'blue');
% hold on
% h_cis = cis_temp.plot('Alpha', 0.1, 'EdgeColor', 'black');
% h_traj = plot3(x_log(1, :), x_log(2, :), x_log(3, :), 'Color', [0 0 0.5]);
% h_dots = scatter3(x_log(1,:), x_log(2,:), x_log(3,:), 30, 'cyan', 'filled');
% 
% title('Traiettoria del sistema termico')
% xlabel('$T_1$ [$^\circ$C]', 'Interpreter', 'latex')
% ylabel('$T_2$ [$^\circ$C]', 'Interpreter', 'latex')
% zlabel('$T_3$ [$^\circ$C]', 'Interpreter', 'latex')
% 
% legend([h_cis, h_npstep, h_traj, h_dots], ...
%     {'CIS', sprintf('%d-step set', N), 'Traiettoria', 'Campioni'}, ...
%     'Interpreter','latex')
% 
% grid on
% view(3)