% Confronto tra le tre implementazioni MPC
% MPC, MPCc e MPC_prova
% Questo script mostra le differenze e confronta le performance

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

%% Simulazione con le tre implementazioni MPC
n_sim = 50; % Ridotto per velocizzare il confronto

% Inizializzazione per MPC classico
htt_mpc = [];
hxx_mpc = [];
u_online_mpc = [];
flag_mpc = zeros(1, n_sim);

% Inizializzazione per MPCc
htt_mpcc = [];
hxx_mpcc = [];
u_online_mpcc = [];
flag_mpcc = zeros(1, n_sim);

% Inizializzazione per MPC_prova
htt_prova = [];
hxx_prova = [];
u_online_prova = [];
flag_prova = zeros(1, n_sim);

% Stato iniziale
x_ini = [284 285 284 0 10 0]';

% Calcolo matrici per MPC e MPCc
[A_cal, A_cal_n, B_cal, B_cal_n, Q_cal, R_cal] = Calligrafica(sys_discreto.A, sys_discreto.B, Q, R, S, Np);

fprintf('=== INIZIO SIMULAZIONE CONFRONTO ===\n');
fprintf('Simulazione con %d iterazioni...\n', n_sim);

for i = 1:n_sim
    
    if i == 1
        x_run = x_ini-x_ref(1:6);
    else
        x_run = hxx_mpc(:, end)-x_ref(1:6);
    end

    % === MPC CLASSICO ===
    [controlAction_mpc, flag_mpc(i)] = MPC(x_run, A_cal, A_cal_n, B_cal, B_cal_n, Q_cal, R_cal, Np, G, g, X_v_lin, U_v_lin);
    
    % === MPCc ===
    [controlAction_mpcc, flag_mpcc(i)] = MPCc(x_run, A_cal, A_cal_n, B_cal, B_cal_n, Q_cal, R_cal, Np, G, g, X_v_lin, U_v_lin);
    
    % === MPC_prova ===
    [controlAction_prova, flag_prova(i)] = MPC_prova(x_run, sys_discreto.A, sys_discreto.B, Q, R, Np, Hx, hx, Hu, hu, G, g, x_ref(1:6), u_ref(1:3));
    
    % Simulazione MPC classico
    tempo = linspace(Ts*(i-1), Ts*i, Ts);
    controlAction_mpc = controlAction_mpc(1:3) + [100; 100; 100];
    u_online_mpc = [u_online_mpc, repmat(controlAction_mpc, 1, Ts)];
    dxdt = @(t,x) CasaEsterno(t, x, k, C, tau, T_ext, k_ext, controlAction_mpc);
    [tt, xx] = ode45(dxdt, tempo, x_run+x_ref);
    htt_mpc = [htt_mpc, tt'];
    hxx_mpc = [hxx_mpc, xx'];
    
    % Simulazione MPCc
    controlAction_mpcc = controlAction_mpcc(1:3) + [100; 100; 100];
    u_online_mpcc = [u_online_mpcc, repmat(controlAction_mpcc, 1, Ts)];
    dxdt = @(t,x) CasaEsterno(t, x, k, C, tau, T_ext, k_ext, controlAction_mpcc);
    [tt, xx] = ode45(dxdt, tempo, x_run+x_ref);
    htt_mpcc = [htt_mpcc, tt'];
    hxx_mpcc = [hxx_mpcc, xx'];
    
    % Simulazione MPC_prova
    controlAction_prova = controlAction_prova(1:3) + [100; 100; 100];
    u_online_prova = [u_online_prova, repmat(controlAction_prova, 1, Ts)];
    dxdt = @(t,x) CasaEsterno(t, x, k, C, tau, T_ext, k_ext, controlAction_prova);
    [tt, xx] = ode45(dxdt, tempo, x_run+x_ref);
    htt_prova = [htt_prova, tt'];
    hxx_prova = [hxx_prova, xx'];
    
    % Progress bar
    if mod(i, 10) == 0
        fprintf('Completato: %d/%d (%.0f%%)\n', i, n_sim, 100*i/n_sim);
    end
end

fprintf('Simulazione completata!\n');

%% Analisi delle performance
tempo_mpc = htt_mpc/60;
tempo_mpcc = htt_mpcc/60;
tempo_prova = htt_prova/60;

% Calcolo statistiche
successi_mpc = sum(flag_mpc);
successi_mpcc = sum(flag_mpcc);
successi_prova = sum(flag_prova);

fprintf('\n=== STATISTICHE CONFRONTO ===\n');
fprintf('MPC classico:  %d/%d successi (%.1f%%)\n', successi_mpc, n_sim, 100*successi_mpc/n_sim);
fprintf('MPCc:          %d/%d successi (%.1f%%)\n', successi_mpcc, n_sim, 100*successi_mpcc/n_sim);
fprintf('MPC_prova:     %d/%d successi (%.1f%%)\n', successi_prova, n_sim, 100*successi_prova/n_sim);

%% Plot di confronto

% Configurazione globale per i plot
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultAxesFontWeight', 'bold');

% Prima figura: Confronto temperature
figure('Name', 'Confronto Temperature', 'NumberTitle', 'off');
sgtitle("Confronto evoluzioni temperature - Tre implementazioni MPC", 'FontSize', 16, 'FontWeight', 'bold');

subplot(2, 1, 1);
plot(tempo_mpc, hxx_mpc(1, :), 'r-', 'LineWidth', 2);  % T1 MPC
hold on;
plot(tempo_mpc, hxx_mpc(2, :), 'g-', 'LineWidth', 2);  % T2 MPC
plot(tempo_mpc, hxx_mpc(3, :), 'b-', 'LineWidth', 2);  % T3 MPC
yline(x_ref(1), 'Color', 'k', 'LineWidth', 2, 'Label', 'Obiettivo');
grid on;
legend(["T1 MPC", "T2 MPC", "T3 MPC", "Obiettivo"], 'Location', 'best');
ylabel("Temperatura $[^{\circ}C]$", "Interpreter", "latex");
title("MPC Classico");
ylim([280, 290]);

subplot(2, 1, 2);
plot(tempo_prova, hxx_prova(1, :), 'r--', 'LineWidth', 2);  % T1 MPC_prova
hold on;
plot(tempo_prova, hxx_prova(2, :), 'g--', 'LineWidth', 2);  % T2 MPC_prova
plot(tempo_prova, hxx_prova(3, :), 'b--', 'LineWidth', 2);  % T3 MPC_prova
yline(x_ref(1), 'Color', 'k', 'LineWidth', 2, 'Label', 'Obiettivo');
grid on;
legend(["T1 MPC_prova", "T2 MPC_prova", "T3 MPC_prova", "Obiettivo"], 'Location', 'best');
ylabel("Temperatura $[^{\circ}C]$", "Interpreter", "latex");
xlabel("Tempo $[min]$", "Interpreter", "latex");
title("MPC_prova");
ylim([280, 290]);

% Seconda figura: Confronto azioni di controllo
figure('Name', 'Confronto Azioni di Controllo', 'NumberTitle', 'off');
sgtitle("Confronto azioni di controllo - Tre implementazioni MPC", 'FontSize', 16, 'FontWeight', 'bold');

subplot(3, 1, 1);
plot(tempo_mpc, u_online_mpc(1, :), 'r-', 'LineWidth', 2);
hold on;
plot(tempo_mpc, u_online_mpc(2, :), 'g-', 'LineWidth', 2);
plot(tempo_mpc, u_online_mpc(3, :), 'b-', 'LineWidth', 2);
grid on;
legend(["Q1", "Q2", "Q3"], 'Location', 'best');
ylabel("Potenza $[W]$");
title("MPC Classico");
ylim([90, 150]);

subplot(3, 1, 2);
plot(tempo_mpcc, u_online_mpcc(1, :), 'r-', 'LineWidth', 2);
hold on;
plot(tempo_mpcc, u_online_mpcc(2, :), 'g-', 'LineWidth', 2);
plot(tempo_mpcc, u_online_mpcc(3, :), 'b-', 'LineWidth', 2);
grid on;
legend(["Q1", "Q2", "Q3"], 'Location', 'best');
ylabel("Potenza $[W]$");
title("MPCc");
ylim([90, 150]);

subplot(3, 1, 3);
plot(tempo_prova, u_online_prova(1, :), 'r-', 'LineWidth', 2);
hold on;
plot(tempo_prova, u_online_prova(2, :), 'g-', 'LineWidth', 2);
plot(tempo_prova, u_online_prova(3, :), 'b-', 'LineWidth', 2);
grid on;
legend(["Q1", "Q2", "Q3"], 'Location', 'best');
ylabel("Potenza $[W]$");
xlabel("Tempo $[min]$", "Interpreter", "latex");
title("MPC_prova");
ylim([90, 150]);

% Terza figura: Confronto flag di successo
figure('Name', 'Confronto Flag di Successo', 'NumberTitle', 'off');
subplot(3, 1, 1);
plot(1:n_sim, flag_mpc, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
title("Flag di successo - MPC Classico");
ylabel("Flag");
ylim([-0.1, 1.1]);
yticks([0, 1]);
yticklabels({'Fallimento', 'Successo'});

subplot(3, 1, 2);
plot(1:n_sim, flag_mpcc, 'g-o', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
title("Flag di successo - MPCc");
ylabel("Flag");
ylim([-0.1, 1.1]);
yticks([0, 1]);
yticklabels({'Fallimento', 'Successo'});

subplot(3, 1, 3);
plot(1:n_sim, flag_prova, 'r-o', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
title("Flag di successo - MPC_prova");
ylabel("Flag");
xlabel("Iterazione");
ylim([-0.1, 1.1]);
yticks([0, 1]);
yticklabels({'Fallimento', 'Successo'});

% Quarta figura: Confronto errori di tracking
figure('Name', 'Confronto Errori di Tracking', 'NumberTitle', 'off');
sgtitle("Confronto errori di tracking - Tre implementazioni MPC", 'FontSize', 16, 'FontWeight', 'bold');

% Calcolo errori
err_temp_mpc = sqrt(mean((hxx_mpc(1:3, :) - x_ref(1:3)).^2, 2));
err_temp_prova = sqrt(mean((hxx_prova(1:3, :) - x_ref(1:3)).^2, 2));

subplot(2, 1, 1);
bar([err_temp_mpc, err_temp_prova]);
grid on;
legend(["MPC Classico", "MPC_prova"], 'Location', 'best');
ylabel("RMSE Temperatura $[^{\circ}C]$", "Interpreter", "latex");
title("Errori di tracking temperature");
xticklabels(["T1", "T2", "T3"]);

% Calcolo errori potenza
err_pot_mpc = sqrt(mean((hxx_mpc(4:6, :) - x_ref(4:6)).^2, 2));
err_pot_prova = sqrt(mean((hxx_prova(4:6, :) - x_ref(4:6)).^2, 2));

subplot(2, 1, 2);
bar([err_pot_mpc, err_pot_prova]);
grid on;
legend(["MPC Classico", "MPC_prova"], 'Location', 'best');
ylabel("RMSE Potenza $[W]$", "Interpreter", "latex");
title("Errori di tracking potenza");
xticklabels(["Q1", "Q2", "Q3"]);

% Reset delle impostazioni globali
set(0, 'DefaultLineLineWidth', 'remove');
set(0, 'DefaultAxesFontSize', 'remove');
set(0, 'DefaultAxesFontWeight', 'remove');

fprintf('\n=== ANALISI COMPLETATA ===\n');
fprintf('I plot mostrano il confronto tra le tre implementazioni MPC.\n');
fprintf('MPC_prova utilizza l''approccio basato su mpc_ingredients.m\n');
fprintf('mentre MPC e MPCc utilizzano l''approccio tradizionale.\n');
