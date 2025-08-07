function x_dot = pendulum(t,x,u,g,l,m,b)
%PENDULUM Modello nello spazio degli stati di un pendolo (angolo e velocità
%angolare)

%   Equazioni di stato
x_dot = zeros(2,1);
x_dot(1) = x(2);
x_dot(2) = -(u/l*cos(x(1)) + g/l*sin(x(1)) + b/(m*l^2)*x(2));
end