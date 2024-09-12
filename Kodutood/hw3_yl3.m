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
%   2.*x.^3 - 11.7.*x.^2 + 17.7.*x - 5 = 0
g1 = @(x) 5.0./(2.*x.^2 - 11.7.*x + 17.7);      % Leiab ainult esimese nullkoha.
dg1dx = @(x) (14.625 - 5.*x)./(8.85 - 5.85.*x + x.^2).^2;
disp(abs(dg1dx(x0)))
disp(abs(dg1dx(x0)) < 1)

g2 = @(x) nthroot(5.0./( 2.0 - 11.7./x + 17.7./x.^2 ), 3); % Leiab ainult esimese nullkoha.
dg2dx = @(x) -(0.57.*(1./(2 + 17.7./x.^2 - 11.7./x)).^(4/3).*(-35.4 + 11.7.*x))./x.^3;
disp(abs(dg2dx(x0)))
disp(abs(dg2dx(x0)) < 1)

g3 = @(x) sqrt(5.0./( 2.0.*x - 11.7 + 17.7./x ));   % Samuti ainult esimene.
dg3dx = @(x) -1.11803.*(2 - 17.7./x.^2).*(1./(-11.7 + 17.7./x + 2.*x)).^(3./2);
disp(abs(dg3dx(x0)))
disp(abs(dg3dx(x0)) < 1)

x_algne = 2.0;  % Võimalikud nullkoha väärtused: x0 = [0.1, 2, 3.5]
tol = 1.0E-6;
max_iter = 100000;

% NOTE: g2 kõige tundlikum. Koondub ainult 0.1 korral. 
% g1 ja g3 rohkem robustsemad ning koonduvad mistahes arvu korral.
[x_uus, count] = him(y, g1, x_algne, tol, max_iter);    % NB!: g1, g2, g3 koonduvad ainult 
fprintf("x_uus = %.4f\n", x_uus)                        % nullkoha väärtusele 0.3651 .
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