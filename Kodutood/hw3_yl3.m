clc
clear

fprintf("\nÜlesanne 3.\n")
% a) graafiline meetod
figure(2)
y = @(x) 2.*x.^3 - 11.7.*x.^2 + 17.7.*x - 5;
fplot(y, [-10, 10])
title("Ül.3 Funktsiooni graafik")
yline(0, 'k-'); % x-telje joon
xlim([-0.5, 5.])
ylim([-10, 10])
grid on

fprintf("\nb) fsolve() lahend\n")
x0 = [0.1, 2, 3.5]; % Punktid, mis visuaalselt umbes läbivad x-telge
x  = fsolve(y,x0);
disp(x)
% Lisan nullkohad eelmisele graafikule
hold on
plot(x, [0,0,0], 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
hold off

fprintf("\nc) Lahend hariliku iteratsioonmeetodiga")
syms x
%   2.*x.^3 - 11.7.*x.^2 + 17.7.*x - 5 = 0
% I variant
g1 = @(x) (5.0 - 2.*x.^3 + 11.7.*x.^2)./17.7;                       % korras. lohakusviga /17.7 jagamisel
% II variant
%g1 = @(x) 5.0./(2.*x.^2 - 11.7.*x + 17.7);                         % Korras. 
dg1dx = matlabFunction(diff((5.0 - 2.*x.^3 + 11.7.*x.^2)./17.7));   % korras
%dg1dx = matlabFunction(diff(5.0./(2.*x.^2 - 11.7.*x + 17.7)));     % Korras.
disp(abs(dg1dx(x0)))
disp(abs(dg1dx(x0)) < 1)

g2 = @(x) sqrt((2.*x.^3 + 17.7.*x - 5.0)./11.7);                    % Korras.
dg2dx = matlabFunction(diff(sqrt((2.*x.^3 + 17.7.*x - 5.0)./11.7)));
disp(abs(dg2dx(x0)))
disp(abs(dg2dx(x0)) < 1)

g3 = @(x) nthroot((5.0 + 11.7.*x.^2 - 17.7.*x)./2.0, 3);            % Korras.
dg3dx = matlabFunction(diff( nthroot((5.0 + 11.7.*x.^2 - 17.7.*x)./2.0, 3) ));
disp(abs(dg3dx(x0)))
disp(abs(dg3dx(x0)) < 1)

x_algne = 3.5;  % Võimalikud nullkoha väärtused: x0 = [0.1, 2, 3.5]
tol = 1.0E-6;
max_iter = 100000;

% NOTE: g1 esimene vaiant peaks koonduma, aga hajub.
[x_uus, count] = him(y, g1, x_algne, tol, max_iter);    % NB!: g1 = @(x) (5.0 - 2.*x.^3 + 11.7.*x.^2)./11.7; ei sobi mingil põhjusel
fprintf("x_uus = %.4f\n", x_uus)                        
fprintf("count = %d\n", count)


% Hariliku iteratsioonimeetodi funktsioon
function [x_uus, count] = him(f, g, x_algne, tol, max_iter)
    % Leiab lahendi hariliku iteratsioonimeetodiga
    % f - algne funktsioon; f(x) = 0.
    % g - iteratsioonifunktsioon
    % x_algne - algväärtus
    % tol - tolerants
    % max_iter - maksimaalne iteratsioonide arv

    % Algväärtused
    x_uus  = x_algne;
    count  = 0;

    % Iteratsiooniprotsess
    while abs(f(x_uus)) >= tol && count < max_iter
        x_uus  = g(x_uus);
        count  = count + 1;
    end

    % Kontroll, kas maksimaalne iteratsioonide arv ületati
    if count >= max_iter
        warning('Maksimaalne iteratsioonide arv ületatud.');
    end
end