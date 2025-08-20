function [H_nsteps,h_nsteps] = controllable_set(Hx,hx,Hu,hu,H_target,h_target,A,B,N)
%CONTROLLABLE_SET Calcolo dell'N-step-controllable set ad un set target
%descritto da H_target * x <= h_target


% inizializzo a partire dal set  più piccolo che posso raggiungere
% L'alogoritmo lavora all'indietro nel tempo, partendo dall'insieme target.
n = size(A,2);
m = size(B,2);

%   Candidato iniziale: dopo 0 passi l'insieme raggiungibile è il target
%   stesso.
H_ii_steps = H_target;
h_ii_steps = h_target;

% si vuole trovare l'insieme di stati x tali che esista un ingresso u
% ammissibile (Hu*u<=hu) per cui lo stato x è ammissibile (Hx*x<= hx) e lo
% staro successivo x_next=A*x + B*u appartiene all'insieme calcolato al
% passo precedente 

for ii=1:N
    %   Computazione in R^(n+m)
    temp = Polyhedron('A',[H_ii_steps*A H_ii_steps*B; ...
        zeros(size(Hu,1),n) Hu],'b',[h_ii_steps; hu]);
    %   Proiezione in R^n
    temp = projection(temp,1:n);

    % rimuovo vincoli ridondanti 
    temp.minHRep(); 
  
    %   Intersezione con X := {x | Hx*x <= hx}
    H_ii_steps = [temp.A; Hx];
    h_ii_steps = [temp.b; hx];
    fprintf("Iterazione %d: numero vincoli = %d\n", ii, size(temp.A,1));
end



H_nsteps = H_ii_steps;
h_nsteps = h_ii_steps;

end

