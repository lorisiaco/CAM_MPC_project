function dxdt = temperature_nonlinear(t, x, u, k, k_ext, T_ext, C, tau)
    % Stato: x = [T1; T2; T3; Q1; Q2; Q3]
    % Ingresso: u = [Q1r; Q2r; Q3r]

    % variabili di stato
    T1 = x(1); T2 = x(2); T3 = x(3);
    Q1 = x(4); Q2 = x(5); Q3 = x(6);

    % ingressi= setpoint termosifoni
    Q1r = u(1); Q2r = u(2); Q3r = u(3);

    % Scambi termici non lineari tra stanze
    k12 = k(1) + 4 / (1 + exp(-0.5 * abs(T1 - T2)));
    k13 = k(1) + 4 / (1 + exp(-0.5 * abs(T1 - T3)));
    k23 = k(2) + 4 / (1 + exp(-0.5 * abs(T2 - T3)));

    % Dinamica delle temperature
    f1 = (Q1 - k12 * (T1 - T2) - k13 * (T1 - T3) - k_ext * (T1 - T_ext)) / C(1);
    f2 = (Q2 + k12 * (T1 - T2) - k23 * (T2 - T3) - k_ext * (T2 - T_ext)) / C(2);
    f3 = (Q3 + k13 * (T1 - T3) + k23 * (T2 - T3) - k_ext * (T3 - T_ext)) / C(3);

    % Dinamica dei termosifoni
    f4 = (Q1r - Q1) / tau(1);
    f5 = (Q2r - Q2) / tau(2);
    f6 = (Q3r - Q3) / tau(3);

    % Derivata dello stato
    dxdt = [f1; f2; f3; f4; f5; f6];
end
