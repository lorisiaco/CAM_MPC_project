clear; clc; close all;
run('inizializzazione.m');

%% Definizione dei vincoli sullo stato e sull'ingresso
% Limiti delle temperature [K]
T_min = 273;  % 0°C
T_max = 323;  % 50°C

% Limiti delle potenze [W]
Q_min = 0;    % 0W
Q_max = 2000; % 2000W

% Matrici dei vincoli Hx*x <= hx e Hu*u <= hu
Hx = [eye(3); -eye(3); zeros(3,3); zeros(3,3)];
hx = [T_max*ones(3,1); -T_min*ones(3,1); Q_max*ones(3,1); -Q_min*ones(3,1)];

Hu = [eye(3); -eye(3)];
hu = [Q_max*ones(3,1); -Q_min*ones(3,1)];

%% Matrici del costo quadratico
Q = 1e4 * eye(6);
R = 1e1 * eye(3);
[~, S] = dlqr(A_d, B_d, Q, R);

%% Calcolo del Control Invariant Set (CIS)
disp('Calcolo del Control Invariant Set...');
[G, g] = CIS(A_d, B_d, zeros(6,1), zeros(3,1), Hx, hx, Hu, hu, Q, R);

disp('Control Invariant Set calcolato con successo');

% Plot del Polyhedron del CIS
try
    CIS_poly = Polyhedron(G, g);
    CIS_poly.minHRep();
    
    % Proiezione sulle prime 3 variabili (temperature)
    CIS_T = projection(CIS_poly, [1,2,3]);
    figure;
    CIS_T.plot();
    title('Proiezione del CIS sulle temperature');
    xlabel('T1 [K]'); ylabel('T2 [K]'); zlabel('T3 [K]');
    hold on;
    x_start_centrato = x_start - x_ref;
    plot3(x_start_centrato(1), x_start_centrato(2), x_start_centrato(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    legend('CIS', 'Stato iniziale centrato');
    hold off;
catch ME
    disp('Errore nel plotting del CIS:');
    disp(ME.message);
end

%% Calcolo dell'N-step controllable set
N = 50; % Numero di passi per il controllable set
disp(['Calcolo del ', num2str(N), '-step controllable set...']);

try
    [Np_steps_H, Np_steps_h, N_steps] = CS(Hx, hx, Hu, hu, G, g, A_d, B_d, x_start - x_ref);
    disp([num2str(N), '-step controllable set calcolato in ', num2str(N_steps), ' passi']);
catch ME
    disp('Errore nel calcolo del controllable set, uso il CIS come fallback:');
    disp(ME.message);
    Np_steps_H = G;
    Np_steps_h = g;
end

%% Preparazione MPC
% Riferimento del sistema linearizzato
x_ref_lin = zeros(6,1);  % Origine, perché traslato x - x_ref
u_ref_lin = zeros(3,1);

% Calcolo matrici per MPC usando Calligrafica
[A_cal, ~, B_cal, ~, Q_cal, R_cal] = Calligrafica(A_d, B_d, Q, R, S, N);

%% Simulazione MPC
T_sim = 80; % Tempo di simulazione

% Matrice in cui salvo la traiettoria degli stati del sistema
x_log = zeros(6, T_sim+1);
% Vettore dove salvo il controllo applicato ad ogni istante
u_log = zeros(3, T_sim);
% Serve per salvare l'esito del quadProg
flags = zeros(1, T_sim);

% Imposto lo stato iniziale del sistema
x_log(:, 1) = x_start;

for tt = 1:T_sim
    % Prendo lo stato attuale linearizzato
    x_current = x_log(:, tt);
    x_lin = x_current - x_ref;
    x_lin_shifted = x_lin - x_ref_lin;

    % QP per MPC
    H_qp = 2 * (B_cal' * Q_cal * B_cal + R_cal);
    f_qp = 2 * x_lin_shifted' * A_cal' * Q_cal * B_cal;

    A_ineq_u = kron(eye(N), Hu);
    b_ineq_u = repmat(hu, N, 1);

    A_ineq_x = [];
    b_ineq_x = [];
    for i = 1:N
        A_state = Hx * A_d^(i-1) * B_d;
        b_state = hx - Hx * A_d^i * x_lin_shifted;
        A_ineq_x = [A_ineq_x; A_state];
        b_ineq_x = [b_ineq_x; b_state];
    end

    A_cis_final = G * A_d^(N-1) * B_d;
    b_cis_final = g - G * A_d^N * x_lin_shifted;

    A_ineq = [A_ineq_u; A_ineq_x; A_cis_final];
    b_ineq = [b_ineq_u; b_ineq_x; b_cis_final];

    % Risolve il problema di ottimizzazione quadratica
    try
        options = optimoptions('quadprog', 'Display', 'off');
        [u_opt, fval, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, [], [], [], [], [], options);
        flags(tt) = exitflag;
        
        if exitflag > 0
            controlAction = u_opt(1:3);
            u_log(:, tt) = u_ref + controlAction;
        else
            warning("Quadprog fallito al passo %d (exitflag = %d), mantengo u precedente.", tt, exitflag);
            if tt > 1
                u_log(:, tt) = u_log(:, tt-1);  % Ripeti l'ultimo controllo valido
            else
                u_log(:, tt) = u_ref;  % Fallback al riferimento
            end
        end
    catch ME
        warning("Errore in quadprog al passo %d: %s", tt, ME.message);
        if tt > 1
            u_log(:, tt) = u_log(:, tt-1);
        else
            u_log(:, tt) = u_ref;
        end
        flags(tt) = -1;
    end

    % Simulazione del sistema
    u_current = u_log(:, tt);
    
    % Dinamica del sistema usando la funzione esistente
    dxdt = @(t, x) CasaEsterno(t, x, k, C, tau, T_ext, k_ext, u_current);
    
    try
        Ts = 60; % Periodo di campionamento [s]
        [~, xx] = ode45(dxdt, [0 Ts], x_current);
        x_log(:, tt+1) = xx(end, :)';
    catch ME
        warning("Errore nella simulazione ODE al passo %d: %s", tt, ME.message);
        % Fallback: stato precedente
        x_log(:, tt+1) = x_current;
    end
end

%% Risultati finali
x_final = x_log(:, end);

fprintf('\n--- Risultati finali della simulazione ---\n');
fprintf('Temperature finali [°C]:\n');
fprintf('  T1 = %.2f\n', x_final(1) - 273);
fprintf('  T2 = %.2f\n', x_final(2) - 273);
fprintf('  T3 = %.2f\n', x_final(3) - 273);

fprintf('\nPotenza termosifoni finali [W]:\n');
fprintf('  Q1 = %.2f\n', x_final(4));
fprintf('  Q2 = %.2f\n', x_final(5));
fprintf('  Q3 = %.2f\n', x_final(6));
fprintf('------------------------------------------\n');

%% Plot dei risultati
% Traslazione del CIS e del set N-step nelle coordinate originali
CIS_shifted = Polyhedron(G, g) + x_ref;
Np_step_set_shifted = Polyhedron(Np_steps_H, Np_steps_h) + x_ref;

% Proiezioni sui soli stati T1, T2, T3
try
    cis_temp = projection(CIS_shifted, [1 2 3]);
    Np_step_temp = projection(Np_step_set_shifted, [1 2 3]);
    
    figure;
    h_npstep = Np_step_temp.plot('Alpha', 0.05, 'LineWidth', 2, 'EdgeColor', 'blue');
    hold on;
    h_cis = cis_temp.plot('Alpha', 0.1, 'EdgeColor', 'black');
    h_traj = plot3(x_log(1, :), x_log(2, :), x_log(3, :), 'Color', [0 0 0.5]);
    h_dots = scatter3(x_log(1,:), x_log(2,:), x_log(3,:), 30, 'cyan', 'filled');
    
    title('Traiettoria del sistema termico');
    xlabel('T_1 [K]');
    ylabel('T_2 [K]');
    zlabel('T_3 [K]');
    
    legend([h_cis, h_npstep, h_traj, h_dots], ...
        {'CIS', sprintf('%d-step set', N), 'Traiettoria', 'Campioni'});
    
    grid on;
    view(3);
catch ME
    disp('Errore nel plotting della traiettoria:');
    disp(ME.message);
end

disp('Script MPC_V_Disugualianza completato con successo!');

 