clear; clc; close all;
run('inizializzazione.m');

%% Parametri MPC
Q = 1e3 * eye(6);
R = 1e1 * eye(3);
[~, S] = dlqr(A_d, B_d, Q, R);
Np = 5;

% Stato iniziale centrato
x_run = x_start - x_ref;

%% Matrici raccolte per MPC
[A_cal, A_cal_n, B_cal, B_cal_n, Q_cal, R_cal] = Calligrafica(A_d, B_d, Q, R, S, Np);

%% Calcolo del Control Invariant Set (CIS)
disp('Calcolo del Control Invariant Set...');
[G, g] = CIS(A_d, B_d, zeros(6,1), zeros(3,1), Hx, hx, Hu, hu, Q, R);

disp('Control Invariant Set calcolato con successo');

% Plot del Polyhedron del CIS (solo se Polyhedron è disponibile)
try
    CIS_poly = Polyhedron(G, g);
    CIS_poly.minHRep();
    % Proiezione sulle prime 3 variabili (temperature)
    CIS_T = CIS_poly.projection(1:3);
    figure;
    CIS_T.plot();
    title('Proiezione del CIS sulle temperature');
    xlabel('T1'); ylabel('T2'); zlabel('T3');
    hold on;
    plot3(x_run(1), x_run(2), x_run(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    legend('CIS', 'Stato iniziale');
    hold off;
catch
    disp('Polyhedron non disponibile: impossibile plottare il CIS.');
end

%% Calcolo dell''N-step controllable set
[Np_steps_H, Np_steps_h, N_steps] = controllable_set(Hx, hx, Hu, hu, G, g, A_d, B_d, x_run);
disp(['N-step controllable set calcolato in ', num2str(N_steps), ' passi']);

%% QP per MPC
H_qp = 2 * (B_cal' * Q_cal * B_cal + R_cal);
f_qp = 2 * x_run' * A_cal' * Q_cal * B_cal;

A_ineq_u = kron(eye(Np), Hu);
b_ineq_u = repmat(hu, Np, 1);

A_ineq_x = [];
b_ineq_x = [];
for i = 1:Np
    A_state = Hx * A_d^(i-1) * B_d;
    b_state = hx - Hx * A_d^i * x_run;
    A_ineq_x = [A_ineq_x; A_state];
    b_ineq_x = [b_ineq_x; b_state];
end

A_cis_final = G * A_d^(Np-1) * B_d;
b_cis_final = g - G * A_d^Np * x_run;

A_ineq = [A_ineq_u; A_ineq_x; A_cis_final];
b_ineq = [b_ineq_u; b_ineq_x; b_cis_final];

try
    options = optimoptions('quadprog', 'Display', 'off');
    [u_opt, fval, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, [], [], [], [], [], options);
    if exitflag > 0
        controlAction = u_opt(1:3);
        flag = 1;
        disp('MPC con vincoli di disuguaglianza risolto con successo');
        disp(['Valore della funzione obiettivo: ', num2str(fval)]);
    else
        controlAction = zeros(3, 1);
        flag = 0;
    end
catch ME
    controlAction = zeros(3, 1);
    flag = 0;
end

if flag == 1
    x_next = A_d * x_run + B_d * controlAction;
    disp('=== RISULTATI MPC ===');
    disp(['Stato attuale (centrato): ', num2str(x_run')]);
    disp(['Controllo ottimale (centrato): ', num2str(controlAction')]);
    disp(['Stato successivo (centrato): ', num2str(x_next')]);
    disp(['Stato attuale (assoluto): ', num2str((x_run + x_ref)')]);
    disp(['Controllo ottimale (assoluto): ', num2str((controlAction + u_ref)')]);
    disp(['Stato successivo (assoluto): ', num2str((x_next + x_ref)')]);
else
    disp('MPC fallito: impossibile trovare una soluzione ammissibile');
end

disp('Script MPC_V_Disugualianza completato'); 