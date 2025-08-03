function [H_nsteps, h_nsteps, Np] = CS(Hx, hx, Hu, hu, G, g, A, B, point)
% CS - Calcolo del n-step controllable set di un sistema lineare
%
% Input:
%   Hx, hx    - Vincoli sullo stato: Hx*x <= hx
%   Hu, hu    - Vincoli sull'ingresso: Hu*u <= hu  
%   G, g      - Regione obiettivo: G*x <= g
%   A, B      - Matrici del sistema discreto
%   point     - Punto di partenza del sistema
%
% Output:
%   H_nsteps  - Matrice vincoli del controllable set
%   h_nsteps  - Vettore vincoli del controllable set
%   Np        - Numero di passi necessari

% Dimensioni del sistema
n = size(A, 2);
m = size(B, 2);

% Inizializzazione del controllable set
H_current = G;
h_current = g;

% Parametri di controllo
max_iterations = 1000;
tolerance = 1e-6;
convergence_check = 10;  % Controlla convergenza ogni 10 iterazioni

fprintf('Inizio calcolo Controllable Set...\n');
tic;

for i = 1:max_iterations
    
    % Calcolo del set ad un passo rispetto a quello precedente
    % x_next = A*x + B*u con vincoli [Hx*x <= hx; Hu*u <= hu]
    
    % Matrice vincoli combinati per il sistema esteso [x; u]
    H_extended = [H_current*A, H_current*B; 
                  zeros(size(Hu, 1), n), Hu];
    h_extended = [h_current; hu];
    
    % Creazione del poliedro esteso
    poly_extended = Polyhedron('A', H_extended, 'b', h_extended);
    
    % Proiezione nello spazio degli stati
    poly_projected = projection(poly_extended, 1:n);
    poly_projected = poly_projected.minHRep();
    
    % Intersezione con i vincoli di stato
    H_new = [poly_projected.A; Hx];
    h_new = [poly_projected.b; hx];
    
    % Creazione del poliedro risultante
    poly_result = Polyhedron('A', H_new, 'b', h_new);
    poly_result = poly_result.minHRep();
    
    % Aggiornamento del controllable set
    H_current = poly_result.A;
    h_current = poly_result.b;
    
    % Controllo se il punto di partenza è contenuto
    if poly_result.contains(point)
        fprintf('Punto di partenza contenuto al passo %d\n', i);
        break;
    end
    
    % Controllo convergenza periodico
    if mod(i, convergence_check) == 0
        fprintf('Passo %d completato\n', i);
        
        % Verifica se il set si è stabilizzato
        if i > convergence_check
            % Calcolo differenza con iterazione precedente
            size_diff = abs(size(H_current, 1) - size(H_prev, 1));
            if size_diff < tolerance
                fprintf('Convergenza raggiunta al passo %d\n', i);
                break;
            end
        end
        H_prev = H_current;
    end
    
    % Controllo dimensione poliedro
    if size(H_current, 1) > 1000
        warning('CS: Poliedro troppo grande, interruzione');
        break;
    end
end

% Tempo di calcolo
elapsed_time = toc;
fprintf('Tempo impiegato: %.2f secondi\n', elapsed_time);

% Risultati finali
Np = i;
H_nsteps = H_current;
h_nsteps = h_current;

% Verifica finale
final_poly = Polyhedron('A', H_nsteps, 'b', h_nsteps);
if final_poly.contains(point)
    fprintf('✓ Controllable set calcolato con successo\n');
    fprintf('  - Numero passi: %d\n', Np);
    fprintf('  - Numero vincoli: %d\n', size(H_nsteps, 1));
else
    warning('CS: Punto di partenza non contenuto nel controllable set');
end

end