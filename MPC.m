function [controlAction, flag] = MPC(x_run, A_cal, A_cal_n, B_cal, B_cal_n, Q_cal, R_cal, Np, G, g, XvincLineari, UvincLineari)
% MPC - Model Predictive Control per sistema di controllo termico
%
% Input:
%   x_run         - Stato attuale del sistema (6x1)
%   A_cal         - Matrice di stato del sistema discreto (6x6)
%   A_cal_n       - Matrice di stato nominale (6x6)
%   B_cal         - Matrice di ingresso del sistema discreto (6x3)
%   B_cal_n       - Matrice di ingresso nominale (6x3)
%   Q_cal         - Matrice di peso per gli stati (6x6)
%   R_cal         - Matrice di peso per gli ingressi (3x3)
%   Np            - Orizzonte di predizione
%   G             - Matrice vincoli combinati
%   g             - Vettore vincoli combinati
%   XvincLineari  - Vincoli lineari sugli stati [max; min]
%   UvincLineari  - Vincoli lineari sugli ingressi [max; min]
%
% Output:
%   controlAction - Azione di controllo ottimale (3x1)
%   flag         - Flag di successo (1=successo, 0=fallimento)

% Dimensioni del sistema
n = size(A_cal, 1);         % Dimensione stato
m = size(B_cal, 2);         % Dimensione ingresso

% Calcolo della matrice Hessiana H
H = 2 * (B_cal' * Q_cal * B_cal + R_cal);

% Impostazione dei vincoli
A_qp = [B_cal;              % Vincolo di massimo dello stato
         -B_cal;             % Vincolo di minimo dello stato
         eye(m);             % Vincolo di massimo dell'ingresso
         -eye(m);            % Vincolo di minimo dell'ingresso
         G * B_cal_n];       % Vincoli combinati

% Definizione vincoli di ingresso e stati centrati
X_max = [];
X_min = [];
U_max = [];
U_min = [];

% Costruzione dei vincoli per l'orizzonte di predizione
for i = 1:Np+1
    X_max = [X_max; XvincLineari(1:2)];
    X_min = [X_min; XvincLineari(3:4)];
    if i ~= Np+1  % Se i è diverso da Np+1
        U_max = [U_max; UvincLineari(1)];
        U_min = [U_min; UvincLineari(2)];
    end
end

% Calcolo del vettore f e b_qp
f = 2 * x_run' * A_cal' * Q_cal * B_cal;

b_qp = [X_max - A_cal * x_run;      % Vincoli stati superiori
         -X_min + A_cal * x_run;     % Vincoli stati inferiori
         U_max;                       % Vincoli ingressi superiori
         -U_min;                      % Vincoli ingressi inferiori
         g - G * A_cal_n * x_run];   % Vincoli combinati
end 