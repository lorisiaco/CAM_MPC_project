function [H_nsteps,h_nsteps] = controllable_set2(Hx,hx,Hu,hu,H_target,h_target,A,B,N)
%CONTROLLABLE_SET Calcolo dell'N-step-controllable set ad un set target
%descritto da H_target * x <= h_target

% 1. Inizializzo a paratire dal set più piccolo che posso raggiungere
% L'algoritmo lavora "all'indietro" nel tempo, partendo dall'insieme target.
n = size(A,2); % Numero di stati
m = size(B,2); % Numero di ingressi

%   Candidato iniziale
H_ii_steps = H_target;
h_ii_steps = h_target;


% 2. Aggiorno iterativamente il set
% Si vuole trovare l'insieme di stati x tali che esista un ingresso 
% u ammissibile (Hu*u <= hu) per cui lo stato x è ammissibile (Hx*x <= hx) 
% e lo stato successivo x_next = A*x + B*u appartiene all'insieme 
% calcolato al passo precedente (X_{k-1}, descritto da H_{k-1}*x_next <= h_{k-1}).
for ii=1:N
    %   Computazione in R^(n+m)
    temp = Polyhedron('A',[H_ii_steps*A H_ii_steps*B; ...
                            zeros(size(Hu,1),n) Hu], ...
                       'b',[h_ii_steps; ...
                            hu]);
    % Questo poliedro temp contiene le coppie (x,u) tali che:
    %   - H_ii_steps*(A*x + B*u) <= h_ii_steps: Lo stato successivo 
    %       A*x + B*u è nell'insieme X_{k-1} (quello definito da 
    %       H_ii_steps, h_ii_steps prima di questa iterazione).
    %   - Hu*u <= hu: L'ingresso u è ammissibile.

    %   Proiezione in R^n
    temp = projection(temp,1:n); % Questa operazione "elimina" u, 
    % mantenendo solo gli stati x per i quali esiste un u che soddisfa 
    % le condizioni di cui sopra. Il risultato è il predecessor set 
    % di un passo Pre(X_{k-1}).
    temp.minHRep(); % Semplifica la rappresentazione del poliedro.
    %   Intersezione con i vincoli di stato X := {x | Hx*x <= hx}
    H_ii_steps = [temp.A; Hx];
    h_ii_steps = [temp.b; hx];
end

H_nsteps = H_ii_steps;
h_nsteps = h_ii_steps;

end

