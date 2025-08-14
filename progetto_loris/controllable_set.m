function [H_nsteps, h_nsteps, Np] = controllable_set(Hx, hx, Hu, hu, H_target, h_target, A, B, point)
% CONTROLLABLE_SET - Calcolo dell'N-step controllable set ad un set target
%
% Input:
%   Hx, hx      - Vincoli sugli stati (Hx*x <= hx)
%   Hu, hu      - Vincoli sugli ingressi (Hu*u <= hu)
%   H_target, h_target - Set target (H_target*x <= h_target)
%   A, B        - Matrici del sistema
%   point       - Punto di partenza per verificare la raggiungibilità
%
% Output:
%   H_nsteps    - Matrice vincoli N-step controllable set
%   h_nsteps    - Vettore vincoli N-step controllable set
%   Np          - Numero di passi necessari

% Inizializzazione
n = size(A, 1);  % Dimensione stato
m = size(B, 2);  % Dimensione ingresso

% Set iniziale
H_curr = H_target;
h_curr = h_target;

% Contatore iterazioni
i = 0;
max_iter = 100;  % Limite massimo iterazioni

while i < max_iter
    i = i + 1;
    
    % Calcolo del set precedente
    H_prev = H_curr;
    h_prev = h_curr;
    
    % Calcolo del nuovo set
    % x(k+1) = A*x(k) + B*u(k)
    % Vincoli: Hx*x(k) <= hx, Hu*u(k) <= hu, H_target*x(k+1) <= h_target
    
    % Eliminazione di u(k) dai vincoli
    % u(k) = B^+ * (x(k+1) - A*x(k)) dove B^+ è la pseudo-inversa di B
    
    % Vincoli combinati
    H_combined = [Hx; Hu * pinv(B)];
    h_combined = [hx; hu];
    
    % Vincoli sul set target
    H_target_ext = H_target * A;
    h_target_ext = h_target;
    
    % Vincoli totali
    H_total = [H_combined; H_target_ext];
    h_total = [h_combined; h_target_ext];
    
    % Creazione del poliedro e minimizzazione
    poly = Polyhedron(H_total, h_total);
    poly = poly.minHRep();
    
    H_curr = poly.A;
    h_curr = poly.b;
    
    % Verifica convergenza
    if isequal(H_curr, H_prev) && isequal(h_curr, h_prev)
        break;
    end
    
    % Verifica se il punto di partenza è raggiungibile
    if all(H_curr * point <= h_curr)
        break;
    end
end

% Output
H_nsteps = H_curr;
h_nsteps = h_curr;
Np = i;

% Verifica finale
if all(H_nsteps * point <= h_nsteps)
    disp(['Punto raggiungibile in ', num2str(Np), ' passi']);
else
    disp('Punto non raggiungibile nel numero massimo di iterazioni');
    Np = inf;
end

end

