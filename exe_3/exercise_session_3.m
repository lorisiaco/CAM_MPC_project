clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
% set(0,'DefaultFigureWindowStyle', 'normal')
set(0,'defaulttextInterpreter','latex')
rng('default');

%%  0. Parametri del pendolo ed equazioni di stato
g = 9.81;
l = 0.5;
m = 0.1;
b = 1e-2;

%   ODE del sistema
dxdt = @(t,x) pendulum(t,x,u,g,l,m,b);
%%  1. Ritratto di fase pendolo inverso senza attrito (autonomo)

b = 0;  % No attrito
u = 0;  % No forzante

%   ODE del sistema autonomo senza attrito
dxdt = @(t,x) pendulum(t,x,u,g,l,m,b);

%   Proviamo a tracciare alcune traiettorie
x0 = [0 0; 0.3 2; pi/2 0]';
figure(1)
for ii = 1:3
    [tt,xx] = ode45(dxdt,[0 10],x0(:,ii));
    plot(atan2(sin(xx(:,1)),cos(xx(:,1))),xx(:,2));
    title('\textbf{Traiettorie senza attrito}')
    xlabel('$\theta$ [rad]');
    ylabel('$\dot{\theta}$ [rad/s]')
    xlim([-pi pi])
    grid on
    hold on
    drawnow
end

% Questa espressione atan2(sin(xx(:,1)),cos(xx(:,1))) ha uno scopo specifico: normalizzare l'angolo del pendolo nell'intervallo [-π, π].
% Per qualsiasi angolo θ, la coppia (sin(θ), cos(θ)) identifica un punto sul cerchio unitario, e atan2(sin(θ), cos(θ)) recupera l'angolo originale, ma nell'intervallo [-π, π].

%   Tracciamo il ritratto di fase
[x1,x2] = meshgrid(-pi:0.3:pi,-5:0.5:5);
x1_dot = x2;
x2_dot = -g/l .* sin(x1);
quiver(x1,x2,x1_dot,x2_dot)
legend({'$x_0 = (0,0)$','$x_0 = (0.3,2)$',...
    '$x_0 = (\pi/2,0)$','Ritratto di fase'},'Interpreter','latex')

%%  2. Ritratto di fase del pendolo con attrito (autonomo)

b = 1e-2;  % No attrito
u = 0;  % No forzante

%   ODE del sistema autonomo senza attrito
dxdt = @(t,x) pendulum(t,x,u,g,l,m,b);

%   Proviamo a tracciare alcune traiettorie
x0 = [0 0; 0.3 2; pi/2 0]';
figure(2)
for ii = 1:3
    [tt,xx] = ode45(dxdt,[0 10],x0(:,ii));
    plot(atan2(sin(xx(:,1)),cos(xx(:,1))),xx(:,2));
    title('\textbf{Traiettorie con attrito}')
    xlabel('$\theta$ [rad]');
    ylabel('$\dot{\theta}$ [rad/s]')
    xlim([-pi pi])
    grid on
    hold on
    drawnow
end

%   Tracciamo il ritratto di fase
[x1,x2] = meshgrid(-pi:0.3:pi,-5:0.5:5);
x1_dot = x2;
x2_dot = -g/l .* sin(x1) - b/(m*l^2)*x2;
quiver(x1,x2,x1_dot,x2_dot)
hold on

legend({'$x_0 = (0,0)$','$x_0 = (0.3,2)$',...
    '$x_0 = (\pi/2,0)$','Ritratto di fase'},'Interpreter','latex')

%%  3. Funzione di Lyapunov per il pendolo inverso con attrito

%   Energia del sistema (modificata)
V = @(x_1,x_2) 1/2 * (m*l^2)*(x_2.^2) + m*g*l*(1-cos(x_1));

%   Verifichiamo sull'ultima simulazione che V decresca
N = size(xx,1);
Vx = zeros(N,1);
for ii = 1:N
    Vx(ii) = V(xx(ii,1),xx(ii,2));
end

figure(3)
plot(Vx)

%   Plottiamo le curve di livello
[x1,x2] = meshgrid(-pi:0.1:pi,-5:0.5:5);
V_contour = V(x1,x2);
figure(2)
contour(x1,x2,V_contour,10,'LineWidth',2);

legend({'$x_0 = (0,0)$','$x_0 = (0.3,2)$',...
    '$x_0 = (\pi/2,0)$','Ritratto di fase','Contour $V(x)$'},'Interpreter','latex')

%%  4. Invariant set per il sistema controllato

%   Sampling time
Ts = 0.1;

%   Sistema linearizzato e discretizzato (Eulero in avanti)
A = [1 Ts; Ts*g/l -Ts*b/(m*l^2)];
B = [0; Ts/l];

%   Vincoli su stato e ingresso
Hx = [eye(2); -eye(2)];
hx = [pi;1;pi;1];
Hu = [2; -2];
hu = 5*ones(2,1);

%   Matrici del costo quadratico
Q = eye(2);
R = 1;

%   Computazione control invariant set
[G,g] = cis(A,B,[0; 0],0,Hx,hx,Hu,hu,Q,R);
CIS = Polyhedron(G,g);
figure(3)
CIS.plot();
title('\textbf{Control invariant set sistema linearizzato}');
xlabel('$\theta$ [rad]');
ylabel('$\dot{\theta}$ [rad/s]');
xlim([-pi,pi]);
ylim([-2,2]);

%%  5. N-step-controllable set dell'invariant set

Np = 10;

[Np_steps_H, Np_steps_h] = controllable_set(Hx,hx,Hu,hu,G,g,A,B,Np);
Np_step_set = Polyhedron('A',Np_steps_H,'b',Np_steps_h);
figure(4)
Np_step_set.plot('Alpha',0);
title('\textbf{CIS e 10-step-controllable set}');
xlabel('$\theta$ [rad]');
ylabel('$\dot{\theta}$ [rad/s]');
hold on;
CIS.plot()
xlim([-pi,pi]);
ylim([-2,2]);

legend({'10-step set','Control-invariant set'})