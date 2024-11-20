clc
clear
% Ülesanne 1
f = @(x) sin(x)./x;
a = pi/4;
b = pi/2;
h = pi/40;

fprintf("Osalõikude arv:\n");
n = (b-a)/h

% integraal vasakpoolse ristkülikvalemiga
% h*f(x0) + h*f(x1)+ ... + h*f(x9)
x_vasak = a:h:b-h;
int_vasak = h*sum(f(x_vasak))

% integraal parempoolse ristkülikvalemiga
% h*f(x1) + h*f(x1)+ ... + h*f(x10)
x_parem = a+h:h:b;
int_parem = h*sum(f(x_parem))

% Ülesanne 2
f = @(x) x.^2.*cos(x.^4);
a = 0;
b = 10;
n = 4000;  % määrame ise
h = (b-a)/n;

fprintf("integraal trapetsvalemiga")
% h/2 * ( f(x0) + 2*f(x1) + 2*f(x2) + 2*f(x3) + f(x4) )
x = a:h:b;
y = f(x);

h/2.0 * ( f(x(1)) + 2.0*f(x(2)) + 2.0*f(x(3)) + 2.0*f(x(4)) + f(x(5)) )
trapz(x,y)

fprintf("Simpsoni valemiga")
h/3.0 * ( f(x(1)) + 4.0*f(x(2)) + 2.0*f(x(3)) + 4.0*f(x(4)) + f(x(5)) )
I_simp = (h/3) * (f(x(1)) + 2*sum(f(x(3:2:end-2))) + 4*sum(f(x(2:2:end-1))) + f(x(end)))
integral(f, a, b)


fprintf("Ülesanne 4.\n");

syms f(x) g(x);
f(x) = x.^4;
int(f(x), x, 0, 1)

g(x) = sin(x)./x;
int(g(x), x, pi/4, pi/2)

a = -1.0;
b = -0.3;
x = -1.0:0.1:-0.3
y = [10 11 11.5 11.6 11.4 11.1 10 8.7]
figure(1)
plot(x, y, "o")
grid on

f = @(xx) interp1(x, y, xx, "spline");
integral(f, a, b)

k2 = polyfit(x, y, 2);
g = @(x) k2(1).*x.^2 + k2(2).*x + k2(3);
hold on
fplot(g, [-1.1, 0])
integral(g, a, b)

fprintf("Ülesanne 5.\n");
x = 0:1:6
F = [0.0, 5.1, 10.2, 16.1, 25.4, 40.2, 61.9]
%       (b
%   A = | F(x) dx
%      a)
figure(2)
plot(x, F, "o") % tundub olevat parabooli kujuga
grid on
coef = polyfit(x, F, 3);
F_func = @(x) coef(1).*x.^3 + coef(2).*x.^2 + coef(3).*x + coef(4);
hold on
fplot(F_func, [0, 6])
a = 0.0;
b = 6.0;
A = integral(F_func, a, b)

fprintf("Ülesanne 6.\n");
x = 1:0.5:5;
y = [0.0, 10.0, 18.0, 25.0, 32.2, 39.0, 43.0, 48.0, 52.3];
a = 1.0;
b = 5.0;
figure(3)
plot(x,y,"o")
grid on

coef = polyfit(x, y, 7);
%f = @(x) coef(1).*x.^6 + coef(2).*x.^5 + coef(3).*x.^4 + coef(4).*x.^3 + coef(5).*x.^2 + coef(6).*x + coef(7);
f = @(x) coef(1).*x.^7 + coef(2).*x.^6 + coef(3).*x.^5 + coef(4).*x.^4 + coef(5).*x.^3 + coef(6).*x.^2 + coef(7).*x.^1 + coef(8);

hold on
fplot(f, [min(x), max(x)])
integral(f, a, b)
