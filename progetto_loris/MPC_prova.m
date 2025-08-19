function [controlAction, flag] = MPC_prova(x_run, A, B, Q, R, Np, Hx, hx, Hu, hu, CIS_H, CIS_h, x_ref, u_ref)
% MPC_prova - Model Predictive Control alternativo basato su mpc_ingredients
%
% Input:
%   x_run         - Stato attuale del sistema (6x1)
%   A             - Matrice di stato del sistema discreto (6x6)
%   B             - Matrice di ingresso del sistema discreto (6x3)
%   Q             - Matrice di peso per gli stati (6x6)
%   R             - Matrice di peso per gli ingressi (3x3)
%   Np            - Orizzonte di predizione
%   Hx            - Matrice vincoli stati
%   hx            - Vettore vincoli stati
%   Hu            - Matrice vincoli ingressi
%   hu            - Vettore vincoli ingressi
%   CIS_H         - Matrice vincoli terminali
%   CIS_h         - Vettore vincoli terminali
%   x_ref         - Stato di riferimento
%   u_ref         - Ingresso di riferimento
%
% Output:
%   controlAction - Azione di controllo ottimale (3x1)
%   flag         - Flag di successo (1=successo, 0=fallimento)

% Dimensioni del problema
n = size(A, 2); % numero degli stati
m = size(B, 2); % numero degli ingressi
n_ter = length(CIS_h); % numero di righe del vincolo terminale

% Matrice per costo terminale (soluzione di Riccati)
[~, P, ~] = dlqr(A, B, Q, R);

% Traslare i vincoli rispetto al riferimento
% Vincoli sullo stato traslato
Hx_shifted = Hx;
hx_shifted = hx - Hx*x_ref;

% Vincoli per l'ingresso
Hu_shifted = Hu;
hu_shifted = hu - Hu*u_ref;

% Peso sugli stati
Q_tilde = kron(eye(Np), Q);
Q_tilde = blkdiag(Q_tilde, P);

% Peso sugli ingressi
R_tilde = kron(eye(Np), R);

% Matrice dipendenza predizioni da stato iniziale
A_cal = zeros(n*(Np+1), n);
for ii = 1:(Np+1)
    if ii == 1
        A_cal((ii-1)*n+1:ii*n, :) = eye(n);
    else
        A_cal((ii-1)*n+1:ii*n, :) = A^(ii-1);
    end
end

% Matrice dipendenza predizioni da ingressi
B_cal = zeros(n*(Np+1), m*Np);
A_cal_times_B = A_cal*B;
for ii = 1:Np
    B_cal(ii*n+1:end, (ii-1)*m+1:ii*m) = A_cal_times_B(1:(Np-ii+1)*n, :);
end

% Matrice hessiana costo quadratico
F = B_cal'*Q_tilde*B_cal + R_tilde;

% Componente lineare costo quadratico
f = B_cal' * Q_tilde * A_cal;

% Vincoli
Hx_tilde = kron(eye(Np+1), Hx_shifted);
hx_tilde = repmat(hx_shifted, [Np+1, 1]);

Hx_tilde = [Hx_tilde; zeros(n_ter, Np*n), CIS_H];
hx_tilde = [hx_tilde; CIS_h];

Hu_tilde = kron(eye(Np), Hu_shifted);
hu_tilde = repmat(hu_shifted, [Np, 1]);

% Set ammissibili di ingresso (inequalities)
A_ineq = [Hx_tilde*B_cal; Hu_tilde];
b_ineq = [hx_tilde; hu_tilde];

% Calcolo del termine costante per i vincoli (dipendente da x_run)
b_ineq_x0 = b_ineq - [Hx_tilde*A_cal; zeros(size(Hu_tilde, 1), n)] * x_run;

% Risoluzione del problema QP
try
    % Usa quadprog se disponibile
    if exist('quadprog', 'file')
        options = optimoptions('quadprog', 'Display', 'off');
        u_opt = quadprog(F, f*x_run, A_ineq, b_ineq_x0, [], [], [], [], [], options);
    else
        % Fallback: minimizzazione diretta con vincoli
        warning('quadprog non disponibile, uso minimizzazione diretta');
        u_opt = -F \ (f*x_run);
    end
    
    % Estrazione della prima azione di controllo
    controlAction = u_opt(1:m);
    flag = 1;
    
catch ME
    % In caso di errore, usa controllo LQR semplificato
    warning('Errore nella risoluzione QP, uso controllo LQR: %s', ME.message);
    K = -dlqr(A, B, Q, R);
    controlAction = K * x_run;
    flag = 0;
end

end
