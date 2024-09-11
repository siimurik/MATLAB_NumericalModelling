clc
clear
% Iseseisev töö nr. 3

disp("Ülesanne 1.")
disp("a)")
% Teist järku tuletise leidmine Symbolic Math Toolbox abil
syms x
y1 = 6*x^3 - x*tan(x);
disp(diff(y1, 2))

disp("b)")
y2 = nthroot(4*x-cos(7*x), 5);
disp(diff(y2, 2))


fprintf("\nÜlesanne 2.\n")
f = @(x) -2.*x.^6 - 1.5.*x.^4 + 10.*x + 2;
x = linspace(-3, 3, 1000);
y = f(x);
disp("b) Graafilise meetodiga reaalarvuliste lahendite leidmine.")
figure(1);
plot(x, y);
grid on;
xlabel('x');
ylabel('f(x)');
title('Funkstiooni graafik');
hold on;
xline(0, 'k-'); % y-telje joon
yline(0, 'k-'); % x-telje joon
xlim([-1.75, 2.25])
ylim([-60, 40])
disp("Graafiku järgi on nullkohad on umbes x1 = -0.20 ja x2 = 1.32.")
hold off

disp("Kontrolliks:")
% f = -2*x.^6 + 0*x^5 - 1.5*x.^4 + 0*x^3 + 0*x^2 + 10*x + 2
coef = [-2 0 -1.5 0 0 10 2];
x_zeros = find_real_roots(coef);

fprintf("\nc) Reaalarvulised lahendid hariliku iteratsioonimeetidiga.\n")
% Funktsiooni tuletamine h.i.m jaoks
%   -2*x^6 - 1.5*x^4 + 10*x + 2 = 0
%   x*(-2*x^5 - 1.5*x^3 + 10) = -2
%   x = -2/(-2*x^5 - 1.5*x^3 + 10)
%   x = g(x)
g = @(x) -2.0/(-2*x^5 - 1.5*x^3 + 10);
x_uus  =  2.0;
x_vana =  1.0;
tol    = 1.0E-6;
count  = 0;
while norm(abs(x_uus-x_vana), 2) >= tol
    x_vana = x_uus;
    x_uus  = g(x_vana);
    count  = count + 1;
end
fprintf('Vastus h.i.m korral: %.4f\n', x_uus)
fprintf('Iteratsioonide arv: %d\n', count)

fprintf("\nÜlesanne 3.\n")
% a) graafiline meetod
figure(2)
y = @(x) 2.*x.^3 - 11.7.*x.^2 + 17.7.*x - 5;
fplot(y, [-10, 10])
yline(0, 'k-'); % x-telje joon
xlim([-0.5, 5.])
ylim([-10, 10])
grid on
hold off

fprintf("\nb) fsolve() lahend\n")
x0 = [0, 2, 3.5]; % Punktid, mis visuaalselt umbes läbivad x-telge
x  = fsolve(y,x0);
disp(x)

fprintf("\nb) Lahend hariliku iteratsioonmeetodiga\n")
%   x*(2*x^2 - 11.7*x + 17.7) = 5
%   x = 5/(2*x^2 - 11.7*x + 17.7)
g = @(x) 5.0/(2*x^2 - 11.7*x + 17.7);
x_algne = 0.0;
tol = 1.0E-6;
[x_uus, count] = harilik_iter(g, x_algne, tol)

fprintf("\nÜlesanne 4.\n")
alpha = 4.0;
beta  = 5.0;
gamma = 1.0;

y = @(x) cos(x) + (-1)^gamma.*0.1.*alpha.*x - 0.001*beta;
figure(3)
fplot(y)
yline(0, 'k-'); % x-telje joon
grid on
hold off

x0 = [1.0]; % Punktid, mis visuaalselt umbes läbivad x-telge
x  = fsolve(y,x0)

fprintf("\nÜlesanne 5.\n")
a = 4.0;
b = 5.0;
c = 1.0;
figure(4)
fplot(@(x) a./10.*x - b./40 - (-1)^c.*sin(x))
grid on

%   a/10*x - b/40 - (-1)^c*sin(x) = 0
%   a/10*x = b/40 + (-1)^c*sin(x)
%   x = 10/a*(b/40 + (-1)^c*sin(x))
g = @(x) 10./a*(b/40. + (-1.0)^c*sin(x));
x_algne = 0.2;
tol = 1.0E-3;
x_uus  = x_algne;
x_vana = x_algne - 1; % Algväärtus, mis erineb x_uus-st
count  = 0;

% Iteratsiooniprotsess
% NOTE: Stuck in infinite loop
while norm(abs(x_uus - x_vana), 2) >= tol
    x_vana = x_uus
    x_uus  = g(x_vana);
    count  = count + 1;
end

function [x_uus, count] = harilik_iter(g, x_algne, tol)
    % Leiab lahendi hariliku iteratsioonimeetodiga
    % g - iteratsioonifunktsioon
    % x_algne - algväärtus
    % tol - tolerants

    % Algväärtused
    x_uus  = x_algne;
    x_vana = x_algne - 1; % Algväärtus, mis erineb x_uus-st
    count  = 0;

    % Iteratsiooniprotsess
    while norm(abs(x_uus - x_vana), 2) >= tol
        x_vana = x_uus;
        x_uus  = g(x_vana);
        count  = count + 1;
    end
end

function x_zeros = find_real_roots(coef)
    % Leiab polünoomi reaalarvulised juured
    % coef - polünoomi koefitsientide vektor
    
    % Leia kõik juured
    r = roots(coef);
    
    % Initsialiseeri tühi vektor reaalarvuliste juurte jaoks
    x_zeros = [];
    
    % Kontrolli iga juure reaalarvulisust
    for i = 1:length(r)
        if imag(r(i)) == 0.0
            fprintf("[%f, %f]\n", real(r(i)), 0.0);
            x_zeros = [x_zeros, real(r(i))];
        end
    end
end