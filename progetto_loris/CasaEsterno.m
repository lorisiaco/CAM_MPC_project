function x_dot = CasaEsterno(t, x, k, C, tau, T_ext, k_ext, u)
% Casa Esterno - Modello non lineare della dinamica termica di una casa con 3 stanze
%
% INPUTS:
%   t      - tempo (richiesto da ODE solver, non usato direttamente)
%   x      - stato del sistema [T1; T2; T3; Q1; Q2; Q3]
%   k      - vettore dei coefficienti termici base [k1; k2; k3]
%   C      - vettore delle capacità termiche [C1; C2; C3]
%   tau    - costanti di tempo dei termosifoni [tau1; tau2; tau3]
%   T_ext  - temperatura esterna (costante)
%   k_ext  - coefficiente di scambio termico con l'esterno
%   u      - input (potenze di riferimento dei termosifoni) [Q1_r; Q2_r; Q3_r]
%
% OUTPUT:
%   x_dot  - derivata dello stato [dT1/dt; dT2/dt; dT3/dt; dQ1/dt; dQ2/dt; dQ3/dt]

    % Inizializzazione
    x_dot = zeros(6,1);         % derivata dello stato
    k_esima = zeros(3,3);       % matrice dei coefficienti k_ij(t)

    % Estrazione variabili di stato
    T = x(1:3);                 % Temperature: T1, T2, T3
    Q = x(4:6);                 % Potenze attuali: Q1, Q2, Q3

    % Calcolo dei coefficienti di scambio termico non lineari k_ij(t)
    for i = 1:3
        for j = 1:3
            if i ~= j
                delta_T = abs(T(i) - T(j));
                k_esima(i,j) = k(i) + 4 / (1 + exp(-0.5 * delta_T));
            end
        end
    end

    % Dinamica termosifoni (Q̇_i)
    for i = 1:3
        x_dot(i+3) = (u(i) - Q(i)) / tau(i);
    end

    % Dinamica temperature (Ṫ_i)
    x_dot(1) = ( Q(1) - k_esima(1,2)*(T(1)-T(2)) - k_esima(1,3)*(T(1)-T(3)) - k_ext*(T(1) - T_ext) ) / C(1);

    x_dot(2) = ( Q(2) + k_esima(1,2)*(T(1)-T(2)) - k_esima(2,3)*(T(2)-T(3)) - k_ext*(T(2) - T_ext) ) / C(2);

    x_dot(3) = ( Q(3) + k_esima(1,3)*(T(1)-T(3)) + k_esima(2,3)*(T(2)-T(3)) - k_ext*(T(3) - T_ext) ) / C(3);
end
