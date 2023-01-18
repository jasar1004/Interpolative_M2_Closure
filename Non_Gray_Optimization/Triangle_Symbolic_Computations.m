clear all
close all
clc

% We aim the find equation for the line passing through (1/3, 1/3), the
% barycentre of the triangle, and any point within the triangle, of
% coordinates (gam1, gam2)
% The equation is of the form: y = a * x  + b

syms gam1 gam2
syms a b

A = [1/3, 1; gam1, 1];
x = [a, b];
c = [1/3; gam2];

% We now aim to solve the system of equations A x = c
x_sol = inv(A)*c;

% Intersection of y = ax + b with y = 1 - x
syms x y a b
eq1 = y == a*x+b;
eq2 = y == 1 - x;

sol = solve(eq1, eq2);

pretty(sol.x)
pretty(sol.y)