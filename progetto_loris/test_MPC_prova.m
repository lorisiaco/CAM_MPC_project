% Test della funzione MPC_prova
% Questo script dimostra come utilizzare la nuova funzione MPC_prova
% basata sul file mpc_ingredients.m

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

%% N-step controllable set
x0_centrato = [284; 285; 284; 0; 10; 0] - x_ref(1:6);
Np = 10;
Np_steps_H = G;
Np_steps_h = g;

disp(" ")
disp("Passi minimi per entrare nel CIS: " + Np);

%% Simulazione a tempo continuo con MPC_prova
htt=[];
hxx = [];
u_online = [];
x_ini = [284 285 284 0 10 0]';
flag = zeros(1, Np);
n_sim = 100;

% Parametri per MPC_prova
% Creiamo le matrici dei vincoli nel formato richiesto
% Hx e hx sono già definiti in inizializzazione
% Hu e hu sono già definiti in inizializzazione

for i = 1:n_sim

    if i == 1
        x_run = x_ini-x_ref(1:6);
    else
        x_run = hxx(:, end)-x_ref(1:6);
    end

    % Chiamata alla nuova funzione MPC_prova
    [controlAction, flag(i)] = MPC_prova(x_run, sys_discreto.A, sys_discreto.B, ...
                                        Q, R, Np, Hx, hx, Hu, hu, G, g, ...
                                        x_ref(1:6), u_ref(1:3));
    
    tempo = linspace(Ts*(i-1), Ts*i, Ts);
    controlAction = controlAction(1:3) + [100; 100; 100];
    u_online = [u_online, repmat(controlAction, 1, Ts)];
    
    dxdt = @(t,x) CasaEsterno(t, x, k, C, tau, T_ext, k_ext, controlAction);
    [tt, xx] = ode45(dxdt, tempo, x_run+x_ref);
    htt = [htt, tt'];
    hxx = [hxx, xx'];

end

%% Plot finale

% Plot diretto senza utilizzare la funzione plotSimulazione
tempo = htt/60; %[min]

% Configurazione globale per i plot
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultAxesFontWeight', 'bold');

% Prima figura: Evoluzioni degli stati
figure('Name', 'Evoluzioni degli stati - MPC_prova', 'NumberTitle', 'off');
sgtitle("Evoluzioni degli stati nel tempo - MPC_prova", 'FontSize', 16, 'FontWeight', 'bold');

% Subplot per le temperature
plot_T = subplot(2, 1, 1);
plot(tempo, hxx(1, :), 'r');  % rosso per T1
hold on;
plot(tempo, hxx(2, :), 'g');  % verde per T2
plot(tempo, hxx(3, :), 'b');  % blu per T3
yline(x_ref(1), 'Color', 'k', 'LineWidth', 2, 'Label', 'Obiettivo');  % nero per obiettivo
grid on;
legend(["T1 - Stanza 1", "T2 - Stanza 2", "T3 - Stanza 3", "Obiettivo"], ...
       'Location', 'best', 'FontSize', 11, 'FontWeight', 'bold');
ylabel("Temperatura $[^{\circ}C]$", "Interpreter", "latex", 'FontSize', 12);
xlabel("Tempo $[min]$", "Interpreter", "latex", 'FontSize', 12);
title("Evoluzione delle temperature nelle stanze - MPC_prova", 'FontSize', 14, 'FontWeight', 'bold');
ylim([280, 290]);
xlim([0, max(tempo)]);
box on;

% Subplot per le potenze termiche
plot_Q = subplot(2, 1, 2);
plot(tempo, hxx(4, :), 'r');  % rosso per Q1
hold on;
plot(tempo, hxx(5, :), 'g');  % verde per Q2
plot(tempo, hxx(6, :), 'b');  % blu per Q3
yline(x_ref(4), 'Color', 'k', 'LineWidth', 2, 'Label', 'Obiettivo');  % nero per obiettivo
grid on;
legend(["Q1 - Termosifone 1", "Q2 - Termosifone 2", "Q3 - Termosifone 3", "Obiettivo"], ...
       'Location', 'best', 'FontSize', 11, 'FontWeight', 'bold');
ylabel("Potenza termica $[W]$", "Interpreter", "latex", 'FontSize', 12);
xlabel("Tempo $[min]$", "Interpreter", "latex", 'FontSize', 12);
title("Evoluzione della potenza termica dei termosifoni - MPC_prova", 'FontSize', 14, 'FontWeight', 'bold');
ylim([0, 120]);
xlim([0, max(tempo)]);
box on;

% Seconda figura: Azioni di controllo
figure('Name', 'Azioni di controllo MPC_prova', 'NumberTitle', 'off');
plot_U = plot(tempo, u_online(1, :), 'r');  % rosso per Q1
hold on;
plot(tempo, u_online(2, :), 'g');  % verde per Q2
plot(tempo, u_online(3, :), 'b');  % blu per Q3
grid on;
title("Azioni di controllo MPC_prova nel tempo", 'FontSize', 16, 'FontWeight', 'bold');
ylabel("Potenza di controllo $[W]$", "Interpreter", "latex", 'FontSize', 12);
xlabel("Tempo $[min]$", "Interpreter", "latex", 'FontSize', 12);
legend(["Q1 - Controllo Stanza 1", "Q2 - Controllo Stanza 2", "Q3 - Controllo Stanza 3"], ...
       'Location', 'best', 'FontSize', 11, 'FontWeight', 'bold');
ylim([90, 150]);
xlim([0, max(tempo)]);
box on;

% Terza figura: Confronto flag di successo
figure('Name', 'Flag di successo MPC_prova', 'NumberTitle', 'off');
plot(1:n_sim, flag, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
title("Flag di successo MPC_prova", 'FontSize', 16, 'FontWeight', 'bold');
ylabel("Flag (1=successo, 0=fallimento)", 'FontSize', 12);
xlabel("Iterazione", 'FontSize', 12);
ylim([-0.1, 1.1]);
yticks([0, 1]);
yticklabels({'Fallimento', 'Successo'});

% Statistiche
successi = sum(flag);
fallimenti = n_sim - successi;
fprintf('\n=== STATISTICHE MPC_prova ===\n');
fprintf('Successi: %d/%d (%.1f%%)\n', successi, n_sim, 100*successi/n_sim);
fprintf('Fallimenti: %d/%d (%.1f%%)\n', fallimenti, n_sim, 100*fallimenti/n_sim);

% Reset delle impostazioni globali
set(0, 'DefaultLineLineWidth', 'remove');
set(0, 'DefaultAxesFontSize', 'remove');
set(0, 'DefaultAxesFontWeight', 'remove');
