clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
% set(0,'DefaultFigureWindowStyle', 'normal')
set(0,'defaulttextInterpreter','latex')
rng('default');

%% Definizione del problema
% Sistema lineare (doppio integratore)
A = [1 1; 0 1];
B = [0.5; 1];

% Definizione dei vincoli sullo stato Hx \leq hx (poliedro nello spazio
% degli stati)
% - 5 <= x_1 <= 5 --> x_1 <= 5 ; -x_1 <= 5
% - 5 <= x_2 <= 5 --> x_2 <= 5 ; -x_2 <= 5
Hx = [eye(2); -eye(2)];
hx = 5*ones(4,1);

% Definizione dei vincoli sugli ingressi Hu \leq hu (poliedro)
% - 1 <= u_1 <= 1 --> u_1 <= 1 ; -u_1 <= 1
Hu = [1; -1];
hu = 1 * ones(2,1);

%% Calcolo e Visualizzazione del Control Invariant Set (CIS)
% Progettazione controllore LQR
Q = eye(2);
R = 1;

% Set point
x_ref = [4.9;0];
u_ref = 0;

[G,g] = cis(A,B,x_ref,u_ref,Hx,hx,Hu,hu,Q,R);

p = Polyhedron('A',G,'b',g); % Questo crea un oggetto Polyhedron 
% (della MPT3 toolbox) che rappresenta il CIS calcolato. 
% La MPT3 toolbox è uno strumento potente per lavorare con poliedri 
% e set theory in MATLAB, essenziale per analisi di sistemi con vincoli.
figure()
hold on
h_cis = p.plot(); % Handle object
h_xref = plot(x_ref(1), x_ref(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'k'); % Ottieni l'handle del punto plottato
xlabel('$x_1$');
ylabel('$x_2$');
legend([h_cis, h_xref], {'Control Invariant Set', 'Riferimento'}, 'Interpreter','latex')

%% Calcolo e Visualizzazione dell'N-Step-Controllable Set
N = 20;

[T,t] = controllable_set(Hx,hx,Hu,hu,G,g,A,B,N);
ctrb_set = Polyhedron('A',T,'b',t);
figure()
h_ctrb=ctrb_set.plot('Alpha',0,'LineWidth',2);
% 'Alpha',0: Imposta la trasparenza del poliedro. 
%  Un valore di 0 significa completamente trasparente, mostrando solo 
%  i bordi. Questo è utile per visualizzare la regione più ampia senza 
%  coprire il CIS interno.
hold on;
h_cis=p.plot();
h_xref = plot(x_ref(1), x_ref(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'k'); % Ottieni l'handle del punto plottato
xlabel('$x_1$');
ylabel('$x_2$');
legend([h_ctrb, h_cis, h_xref], ...
       {sprintf('%d-step set', N), 'Control-invariant set', 'Riferimento'}, ... % Usa sprintf per includere il valore di N
       'Interpreter', 'latex');

% La proprietà fondamentale è che il Control-invariant set (CIS) è sempre contenuto all'interno dell'N-step-controllable set verso se stesso

