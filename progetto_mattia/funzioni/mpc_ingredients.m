function mpc = mpc_ingredients(A,B,Hx,hx,Hu,hu,CIS_H,CIS_h,x_ref,u_ref,Q,R,Np)

% Dimensioni del problema
n = size(A, 2); % numero degli stati
m = size(B, 2); % numero degli ingressi
n_ter = length(CIS_h); % numero di righe del vincolo terminale

% Matrice per costo terminale
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

%   Matrice dipendenza predizioni da stato iniziale
A_cal = zeros(n*(Np+1), n);
for ii = 1:(Np+1)
    if ii == 1
        A_cal((ii-1)*n+1:ii*n, :) = eye(n);
    else
        A_cal((ii-1)*n+1:ii*n, :) = A^(ii-1);
    end
end

%   Matrice dipendenza predizioni da ingressi
B_cal = zeros(n*(Np+1), m*Np);
A_cal_times_B = A_cal*B;
for ii = 1:Np
    B_cal(ii*n+1:end, (ii-1)*m+1:ii*m) = A_cal_times_B(1:(Np-ii+1)*n, :);
end

%   Matrice hessiana costo quadratico
F = B_cal'*Q_tilde*B_cal + R_tilde;

%   Componente lineare costo quadratico
f = B_cal' * Q_tilde * A_cal;

% Vincoli
Hx_tilde = kron(eye(Np+1), Hx_shifted);
hx_tilde = repmat(hx_shifted, [Np+1, 1]);

Hx_tilde = [Hx_tilde; zeros(n_ter, Np*n), CIS_H];
hx_tilde = [hx_tilde; CIS_h];

Hu_tilde = kron(eye(Np), Hu_shifted);
hu_tilde = repmat(hu, [Np, 1]);

% Set ammissibili di ingresso (inequalities)
A_ineq = [Hx_tilde*B_cal; Hu_tilde];
b_ineq = [hx_tilde; hu_tilde];

%   Creazione struttura mpc
mpc.F = F;
mpc.f_base = f;
mpc.A_ineq = A_ineq;
mpc.b_ineq_base = b_ineq;
mpc.Np = Np;
mpc.b_ineq_x0_factor = [Hx_tilde*A_cal; zeros(2*m*Np, n)];
end

