

clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
set(0,'defaulttextInterpreter','latex')
rng('default');

%% Invariant set per il sistema controllato

addpath('funzioni');
main;

% Matrici del costo quadratico
Q = 1.e3*eye(6);
R = 1e1*eye(3);                
% Control invariant set CIS_H*x <= CIS_h
[CIS_H, CIS_h] = cis(A_d, B_d, zeros(6,1), zeros(3,1), Hx, hx, Hu, hu, Q, R); 

%% Plot del CIS

CIS_G = Polyhedron(CIS_H, CIS_h);

% minimal H representation (elimina i vincoli ridondanti, ossia quelli che
% non  influenzano il set. Il risultato è equivalente ma con il minor
% numero di disuguaglianze lineari.
CIS_G = CIS_G.minHRep();
% crea una proiezione del poliedro CIS_G sul piano formato dalle variabili 
% 1 e 3 dello stato. Quindi, se lo stato x appartiene a R^6, sto
% proiettando il poliedro sui soli assi x1 e x3, ad esempio.
% RICORDA CHE: x = [T1 T2 T3 Q1 Q2 Q3]'
CIS_G_H13 = projection(CIS_G, [1,2,3]); %cubo T1 T2 T3
CIS_G_H24 = projection(CIS_G, [4,5,6]); % cubo Q1 Q2 Q3

dim = CIS_G.Dim;
disp(['La dimensione del poliedro è: ', num2str(dim)]);

% Plot delle proiezioni
figure

subplot(1 , 2 , 1)
CIS_G_H13.plot();

subplot(1 , 2 , 2)
CIS_G_H24.plot();

