clear; 
clc; 
close all

%% %% Tempo di Simulazione
Ts = 60; % [secondi] - tempo di simulazione
%% Richiamo file inizializzazione per le variabili
inizializzazione

%% Definizione Q,R ed S
Q = 1.e4 * eye(6);
R = 1e1 * eye(3);

[~, S] = dlqr(sys_discreto.A, sys_discreto.B, Q, R); % S da Riccati

%% Vincoli di uguaglianza
%x = 0 (vincolo di uguaglianza)
G = [eye(6);
     -eye(6)];
g = [zeros(12,1)];

setIniziale = Polyhedron(G, g);
disp("Il set iniziale contiene l'origine? " + setIniziale.contains(zeros(6,1)))

%% Esistenza del Set Terminale di Uguaglianza
% Per i vincoli di uguaglianza, il set terminale è semplicemente l'origine
% Non è necessario calcolare un CIS come nel caso di disuguaglianza

% I plot del Set Terminale di Uguaglianza sono rimossi per semplicità
% Il set terminale è il punto origine (0,0,0) per i vincoli di uguaglianza

%% Calcolo del CIS per confronto (opzionale)
% Anche se non necessario per i vincoli di uguaglianza, calcoliamo il CIS
% per confronto con il caso di disuguaglianza
try
    [G_CIS, g_CIS] = CIS(sys_discreto.A, sys_discreto.B, zeros(6,1), zeros(3,1), Hx, hx, Hu, hu, Q, R);
    
    %% Plot del CIS per confronto
    CIS_G = Polyhedron(G_CIS, g_CIS);
    CIS_G = minHRep(CIS_G);
    CIS_G_T = projection(CIS_G, 1:3);
    CIS_G_Q = projection(CIS_G, 4:6);
    
    figure
    subplot(1, 2, 1)
    CIS_G_T.plot('Color', [0.2, 0.6, 0.8]); % Colore blu elegante
    title("CIS per Confronto - Temperature")
    limitiTemp = [X_v_lin(7), X_v_lin(1)];
    xlim(limitiTemp); ylim(limitiTemp); zlim(limitiTemp);
    xlabel("T1 $[^{\circ}C]$", "Interpreter", "latex")
    ylabel("T2 $[^{\circ}C]$", "Interpreter", "latex")
    zlabel("T3 $[^{\circ}C]$", "Interpreter", "latex")
    
    subplot(1, 2, 2)
    CIS_G_Q.plot('Color', [0.2, 0.6, 0.8]); % Colore blu elegante
    title("CIS per Confronto - Potenze Termiche")
    limitiQ = [X_v_lin(10), X_v_lin(4)];
    xlim(limitiQ); ylim(limitiQ); zlim(limitiQ);
    xlabel("Q1 $[W]$", "Interpreter", "latex")
    ylabel("Q2 $[W]$", "Interpreter", "latex")
    zlabel("Q3 $[W]$", "Interpreter", "latex")
    
    disp("CIS calcolato con successo per confronto")
catch ME
    disp("Impossibile calcolare il CIS per confronto: " + ME.message)
end

%% N-step controllable set 
x0_centrato = [284; 285; 284; 0; 10; 0] - x_ref(1:6);
[Np_steps_H, Np_steps_h] = controllable_set(Hx, hx, Hu, hu, G, g, A_d, B_d, 50);

% Impostiamo manualmente il numero di passi
Np = 20; % Per vincoli di uguaglianza, Np deve essere almeno 20
disp(" ")
disp("Numero di passi per MPC: " + Np);

%% Verifica fattibilità dal punto di partenza
% Per i vincoli di uguaglianza, la verifica è più semplice
x0_centrato = [284; 285; 284; 0; 10; 0] - x_ref(1:6);
disp("Punto di partenza centrato: " + mat2str(x0_centrato(1:3), 1) + " (temperature)");
disp("Punto di partenza centrato: " + mat2str(x0_centrato(4:6), 1) + " (potenze)");

%% Calcolo N-Step CS per confronto (opzionale)
% Anche se non necessario per i vincoli di uguaglianza, calcoliamo il N-Step CS
% per confronto con il caso di disuguaglianza
try
    [Np_steps_H_comp, Np_steps_h_comp] = controllable_set(Hx, hx, Hu, hu, G_CIS, g_CIS, A_d, B_d, 50);
    
    %% Plot del N-Step CS per confronto
    Np_step_comp = Polyhedron(Np_steps_H_comp, Np_steps_h_comp);
    Np_step_comp = Np_step_comp.minHRep();
    Np_steps_T_comp = projection(Np_step_comp, 1:3);
    Np_steps_Q_comp = projection(Np_step_comp, 4:6);
    
    figure
    subplot(1, 2, 1)
    h1 = Np_steps_T_comp.plot('Color', [0.2, 0.6, 0.8]); % Colore blu elegante
    alpha(h1, 0.3); % Imposta trasparenza
    title("N-Step CS per Confronto - Temperature")
    xlim(limitiTemp); ylim(limitiTemp); zlim(limitiTemp);
    hold on
    % Punto di partenza
    plot3(x0_centrato(1), x0_centrato(2), x0_centrato(3), '.', 'MarkerSize', 50, 'Color', 'r', 'DisplayName', 'Start Point')
    xlabel("T1 $[^{\circ}C]$", "Interpreter", "latex")
    ylabel("T2 $[^{\circ}C]$", "Interpreter", "latex")
    zlabel("T3 $[^{\circ}C]$", "Interpreter", "latex")
    legend('N-Step CS', 'Start Point', 'Location', 'northeast', 'FontSize', 10)
    hold off
    
    subplot(1, 2, 2)
    h2 = Np_steps_Q_comp.plot('Color', [0.2, 0.6, 0.8]); % Colore blu elegante
    alpha(h2, 0.3); % Imposta trasparenza
    title("N-Step CS per Confronto - Potenze Termiche")
    xlim(limitiQ); ylim(limitiQ); zlim(limitiQ);
    hold on
    % Punto di partenza
    plot3(x0_centrato(4), x0_centrato(5), x0_centrato(6), '.', 'MarkerSize', 50, 'Color', 'r', 'DisplayName', 'Start Point')
    xlabel("Q1 $[W]$", "Interpreter", "latex")
    ylabel("Q2 $[W]$", "Interpreter", "latex")
    zlabel("Q3 $[W]$", "Interpreter", "latex")
    legend('N-Step CS', 'Start Point', 'Location', 'northeast', 'FontSize', 10)
    hold off
    
    disp("N-Step CS calcolato con successo per confronto")
catch ME
    disp("Impossibile calcolare il N-Step CS per confronto: " + ME.message)
end

%% Simulazione a tempo continuo con MPC
htt=[];
hxx = [];
u_online = [];
flag = [];
x_ini = [284 285 284 0 10 0]';
n_sim = 100;

% Calcolo delle matrici nella forma raccolta
[A_cal, A_cal_n, B_cal, B_cal_n, Q_cal, R_cal] = Calligrafica(sys_discreto.A, sys_discreto.B, Q, R, S, Np);

for i = 1:n_sim

    if i == 1
        x_run = x_ini-x_ref(1:6);
    else
        x_run = hxx(:, end)-x_ref(1:6);
    end

    % Utilizzo della funzione MPCu per vincoli di uguaglianza
    [controlAction, flag(i)] = MPCu(x_run, A_cal, A_cal_n, B_cal, B_cal_n, Q_cal, R_cal, Np, G, g, X_v_lin, U_v_lin);
    
    tempo = linspace(Ts*(i-1), Ts*i, Ts);
    controlAction = controlAction(1:3) + [100; 100; 100];
    u_online = [u_online, repmat(controlAction, 1, Ts)];
    
    % Simulazione del sistema dinamico
    dxdt = @(t,x) CasaEsterno(t, x, k, C, tau, T_ext, k_ext, controlAction);
    [tt, xx] = ode45(dxdt, tempo, x_run+x_ref);
    htt = [htt, tt'];
    hxx = [hxx, xx'];

end

%% Plot finale

tempo = htt/60; %[min]

% Configurazione globale per i plot
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultAxesFontWeight', 'bold');

% Prima figura: Evoluzioni degli stati
figure('Name', 'Evoluzioni degli stati nel tempo', 'NumberTitle', 'off');
sgtitle("Evoluzioni degli stati nel tempo - MPC con Vincoli di Uguaglianza", 'FontSize', 16, 'FontWeight', 'bold');

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
title("Azioni di controllo MPC nel tempo - Vincoli di Uguaglianza", 'FontSize', 16, 'FontWeight', 'bold');
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

% Display dei risultati finali
disp(" ")
disp("=== RISULTATI SIMULAZIONE MPC CON VINCOLI DI UGUAGLIANZA ===");
disp("Numero di iterazioni: " + n_sim);
disp("Iterazioni convergenti: " + sum(flag == 1));
disp("Iterazioni non convergenti: " + sum(flag == -1));
disp("Tasso di convergenza: " + round(100*sum(flag == 1)/n_sim, 2) + "%");
disp("Stato finale centrato: " + mat2str(round(hxx(:, end) - x_ref(1:6), 3), 1));
disp("Errore finale: " + round(norm(hxx(:, end) - x_ref(1:6)), 4));
