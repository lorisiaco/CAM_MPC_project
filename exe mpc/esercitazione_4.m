clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
% set(0,'DefaultFigureWindowStyle', 'normal')
set(0,'defaulttextInterpreter','latex')
rng('default');

%% 0. Parametri del pendolo
g = 9.81;
l = 0.5;
m = 0.1;
b = 1e-2;

%% 1. Invariant set per il sistema controllato
% Tempo di campionamento
Ts = 0.1;

% Sistema linearizzato e discretizzato
A = [1 Ts;
    Ts*g/l 1-Ts*b/(m*l^2)];
B = [0;
    Ts/l];

% Vincoli su stao e ingresso (nelle coordinate del sistema linearizzato)
% -pi/2 <= x1 <= pi/2 --> x1 <= pi/2        -x1 <= pi/2
% -2 <= x2 <= 2       --> x2 <= 2           -x2 <= 2
% -5 <= u <= 5        --> u <= 5            -u <= 5
Hx = [eye(2); -eye(2)];
hx = [pi/2; 2; pi/2; 2];
Hu = [1; -1];
hu = [5; 5]; % 5*ones(2,1);

% Matrici del costo quadratico
Q = eye(2);
R = 1;

% Control invariant set CIS_H*x <= CIS_h
[CIS_H, CIS_h] = cis(A, B, [0; 0], 0, Hx, hx, Hu, hu, Q, R);
CIS = Polyhedron(CIS_H, CIS_h);
figure(1)
h_cis = CIS.plot();
title('\textbf{Control Invariant Set del sistema linearizzato}')
xlabel('$\theta$ [rad]')
ylabel('$\dot{\theta}$ [rad/s]')
xlim([-pi pi])
ylim([-2 2])
legend(h_cis, {'CIS'}, 'Interpreter','latex')

%%  2. N-step-controllable set dell'invariant set
% Orizzonte di predizione
N = 5;

[Np_steps_H, Np_steps_h] = controllable_set(Hx, hx, Hu, hu, CIS_H, CIS_h, A, B, N);
Np_steps_set = Polyhedron('A', Np_steps_H, 'b', Np_steps_h);
figure(2)
h_cis = CIS.plot();
hold on
h_nsteps = Np_steps_set.plot('Alpha', 0, 'LineWidth', 2);
title(sprintf('\textbf{CIS e %d-step-controllable set del sistema linearizzato}', N))
xlabel('$\theta$ [rad]')
ylabel('$\dot{\theta}$ [rad/s]')
xlim([-pi pi])
ylim([-2 2])
legend([h_cis, h_nsteps], ...
    {'CIS', sprintf('%d-step set', N)},...
    'Interpreter','latex')

%% 3. MPC e simulazione
% Tempo di simulazione
T_sim = 60;

% Riferimento
x_ref = [pi; 0];
u_ref = 0;

% Riferimento del sistema linearizzato
x_ref_lin = [0; 0];
u_ref_lin = 0;

mpc = mpc_ingredients(A,B,Hx,hx,Hu,hu,CIS_H,CIS_h,x_ref_lin,u_ref_lin,Q,R,N);

% Log stati e ingressi del sistema
x_log = zeros(2, T_sim+1);
u_log = zeros(1, T_sim);
flags = zeros(1, T_sim);

x_log(:, 1) = [pi-0.4; -0.2];

for tt = 1:T_sim
    % Stato del sistema linearizzato
    x_current = x_log(:, tt);
    x_lin = x_current - [pi; 0];
    x_lin_shifted = x_lin - x_ref_lin;

    % Impostare i vincoli MPC in base alla condizione iniziale
    f = mpc.f_base * x_lin_shifted;
    b_ineq = mpc.b_ineq_base - mpc.b_ineq_x0_factor*x_lin_shifted;

    % Risoluzione del problema di ottimizzazione
    [delta_u_seq, ~, exitflag] = quadprog(mpc.F, f, mpc.A_ineq, b_ineq);

    flags(tt) = exitflag;
    delta_u_seq_first = delta_u_seq(1);

    % Avanzamento/Risposta del sistema
    u_log(tt) = u_ref + delta_u_seq_first;

    dxdt = @(t,x) pendulum(t, x, u_log(tt), g, l, m, b);

    [~, xx] = ode45(dxdt, [0 Ts], x_current);

    x_log(:, tt+1) = xx(end, :)';
end

%% 4. Plot dei risultati
%   Traslazione del CIS e dell'N-step set nelle coordinate originali
CIS_shifted = CIS + x_ref;
Np_step_set_shifted = Np_steps_set + x_ref;

figure(3)
h_npstep_shited = Np_step_set_shifted.plot('Alpha', 0, 'LineWidth', 2);
hold on
h_cis_shifted = CIS_shifted.plot();
h_traj = plot(x_log(1, :), x_log(2, :), 'Color',[0 0 0.5]);
h_traj_dots = scatter(x_log(1,:),x_log(2,:),'cyan');
title('Traiettoria del sistema')
xlabel('$\theta$ [rad]')
ylabel('$\dot{\theta}$ [rad/s]')
xlim([0 2*pi])
ylim([-2 2])

legend([h_cis_shifted, h_npstep_shited, h_traj, h_traj_dots], ...
    {'CIS', sprintf('%d-step set', N), 'Traiettoria', 'Sample'},...
    'Interpreter','latex')






