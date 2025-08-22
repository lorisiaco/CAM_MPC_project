function [G,g] = cis(A,B,x_ref,u_ref,Fx,fx,Fu,fu,Q,R)
%CIS Calcolo del control invariant set (CIS) di un sistema lineare
%   Questo metodo assume che un controllore LQR venga utilizzato
%   all'interno del CIS
%   Input:
%       - A, B: matrici del sistema
%       - x_ref: equilibrio attorno al quale costruire il CIS
%       - Fx*x<=fx: vincoli sullo stato
%       - Fu*u<=fu: vincoli sull'ingresso
%       - Q,R: matrici per LQR

%   1. Controllore LQR --> Viene calcolato il guadagno K del controllore LQR. 
%   La legge di controllo che si assume attiva all'interno del CIS è u(k) = K*(x(k) - x_ref) + u_ref.
K = -dlqr(A,B,Q,R);

%   2. Matrice A del sistema controllato con LQR
A_lqr = A+B*K;

%   3. Vincoli sullo stato e sull'ingresso - combinati (F*x <= f)
F = [Fx; Fu*K];
f = [fx; fu + Fu*(K*x_ref - u_ref)];

%   4. Calcolo del CIS (G*x<=g)
CIS_poly_prev = Polyhedron();
CIS_poly_curr = Polyhedron(F,f);

while CIS_poly_prev.isEmptySet || CIS_poly_prev ~= CIS_poly_curr
    %   CIS_poly_prev ~= CIS_poly_curr: L'iterazione continua finché il candidato CIS non converge.

    %   Memorizza vecchio candidato
    CIS_poly_prev = CIS_poly_curr;
    
    %   Calcola nuovo candidato (G_hat * x <= g_hat)
    %   Si calcola un nuovo candidato CIS_poly_curr. Questo nuovo insieme deve contenere stati x tali che
    %       - x stesso soddisfa i vincoli originali (F*x <= f).
    %       - Lo stato successivo x_next = A_lqr*x + B*(u_ref - K*x_ref) 
    %         (che chiameremo A_cl*x + x_offset) deve appartenere al precedente 
    %         candidato CIS_poly_prev. Questo si traduce nella condizione 
    %         CIS_poly_prev.A * (A_lqr*x + B*(u_ref - K*x_ref)) <= CIS_poly_prev.b. 
    %         Riscrivendo: (CIS_poly_prev.A * A_lqr)*x <= CIS_poly_prev.b - CIS_poly_prev.A * B*(u_ref - K*x_ref).

    G_hat = [CIS_poly_curr.A * A_lqr; F];
    g_hat = [CIS_poly_prev.b + CIS_poly_curr.A*B*(K*x_ref - u_ref); f];
    CIS_poly_curr = Polyhedron(G_hat,g_hat);
    
end

%   5. Disequazioni che descrivono il CIS
% Le matrici G e g descrivono il CIS calcolato come un poliedro G*x <= g.
G = CIS_poly_curr.A;
g = CIS_poly_curr.b;