function [controlAction, flag] = MPCu(x_attuale, A_cal, A_cal_n, B_cal, B_cal_n, Q_cal, R_cal, Window, G, g, Vinc_X, Vinc_U)
% MPCu - Controllo MPC con vincoli di uguaglianza terminali
% 
% Questa funzione implementa un controllore MPC che risolve un problema di
% ottimizzazione quadratica con vincoli di uguaglianza sul set terminale.
% Il controllo è progettato per garantire che lo stato finale soddisfi
% i vincoli G*x_finale = g.
%
% Input:
%   x_attuale: stato corrente del sistema (vettore colonna)
%   A_cal, B_cal: matrici del sistema nella forma raccolta
%   A_cal_n, B_cal_n: matrici per l'orizzonte finale
%   Q_cal, R_cal: matrici di costo per stati e ingressi
%   Window: orizzonte di predizione
%   G, g: vincoli di uguaglianza terminali (G*x = g)
%   Vinc_X, Vinc_U: vincoli su stati e ingressi
%
% Output:
%   controlAction: sequenza ottimale di controlli
%   flag: flag di convergenza dell'ottimizzazione

% Costruzione della matrice Hessiana per il problema QP
% La funzione obiettivo è: J = (1/2)*u'*H*u + f'*u
H = 2 * (B_cal' * Q_cal * B_cal + R_cal);

% Definizione dei vincoli di disuguaglianza
% Vincoli sugli stati: B_cal*u <= X_max - A_cal*x e -B_cal*u <= -X_min + A_cal*x
% Vincoli sugli ingressi: u <= U_max e -u <= -U_min
A_qp = [B_cal;           % vincoli superiori sugli stati
        -B_cal;          % vincoli inferiori sugli stati  
        eye(width(B_cal)); % vincoli superiori sugli ingressi
        -eye(width(B_cal))]; % vincoli inferiori sugli ingressi

% Vincoli di uguaglianza terminali: G * x_finale = g
% Dove x_finale = A_cal_n * x_attuale + B_cal_n * u
A_eq = G * B_cal_n;

% Costruzione dei vettori dei limiti per tutti i passi dell'orizzonte
X_max = [];
X_min = [];
U_max = [];
U_min = [];

for step = 1:Window+1
   % Vincoli sugli stati per ogni passo dell'orizzonte
   X_max = [X_max; Vinc_X(1:6)];
   X_min = [X_min; Vinc_X(7:end)];
   
   % Vincoli sugli ingressi (non per l'ultimo passo)
   if step ~= Window+1
        U_max = [U_max; Vinc_U(1:3)];
        U_min = [U_min; Vinc_U(4:end)];
   end
end

% Vettore gradiente della funzione obiettivo
f = 2 * x_attuale' * A_cal' * Q_cal * B_cal;

% Vettore dei termini noti per i vincoli di disuguaglianza
b_qp = [X_max - A_cal * x_attuale;    % vincoli superiori sugli stati
         -X_min + A_cal * x_attuale;   % vincoli inferiori sugli stati
         U_max;                         % vincoli superiori sugli ingressi
         -U_min];                       % vincoli inferiori sugli ingressi

% Termine noto per i vincoli di uguaglianza terminali
b_eq = g - G * A_cal_n * x_attuale;

% Risoluzione del problema di programmazione quadratica
% Il solver quadprog trova il minimo di (1/2)*u'*H*u + f'*u
% soggetto ai vincoli A_qp*u <= b_qp e A_eq*u = b_eq
[controlAction, ~, flag] = quadprog(H, f, A_qp, b_qp, A_eq, b_eq);

end
