clc
clear

fprintf("Ülesanne 1.\n");
x = [1.0, 1.1, 1.2,  1.3, 1.4, 1.5, 1.6, 1.7,  1.8,  1.9,  2.0];
y = [4.0, 7.0, 2.1, 15.0, 3.0, 5.2, 6.0, 8.5, 14.0, 15.1, 17.8];
h = 0.1; % sammu pikkus

figure(1)
plot(x, y, "o")
grid on
title("Ülesanne 1")


dy = zeros(1, length(x));
% Tuletis esimeses sõlmes
dy(1) = (y(2) - y(1))/h;
% Tuletised sisesõlmedes
for i = 2:(length(x)-1)
    dy(i) = (y(i+1) - y(i-1))/(2*h);
end
% Tuletis viimases sõlmes
dy(length(x)) = ( y(length(x)) - y(length(x)-1) )/h;
disp(dy)

fprintf("\nÜlesanne 2.\n");
x = [7.22, 7.26, 7.30, 7.34];
y = [2.35, 2.29, 2.25, 2.28];
h = 0.04;

figure(2)
plot(x, y, "o")
grid on
title("Ülesanne 2")

% Tuletis 2. sõlmes
(y(3) - 2*y(2) + y(1))/(h*h)
% Tuletis 3. sõlmes
(y(4) - 2*y(3) + y(2))/(h*h)

fprintf("\nÜlesanne 3.\n");
x  = [-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3];
y  = [10.0 11.0 11.5 11.6 11.4 11.1 10.0  8.7];
h = x(2) - x(1);

dy = zeros(1, length(y));
dy(1) = (y(2) - y(1))/h;
for i = 2:(length(x)-1)
    dy(i) = (y(i+1) - y(i-1))/(2*h);
end
dy(length(x)) = ( y(length(x)) - y(length(x)-1) )/h;
dy

ddy = zeros(1, length(y));
for i = 2:(length(x)-1)
    ddy(i) = (y(i+1) - 2*y(i) + y(i-1))/(h*h);
end
ddy

fprintf("b) lähendame esialgselt funkstiooni kuupsplainiga S^(3,2)\n")
xx = min(x): 0.001: max(x);
S3_2 = interp1(x, y, xx, "spline");
hh = 0.001;

figure(3)
plot(x, y, "o")
grid on
hold on
plot(xx, S3_2, 'm')
plot(x, dy, "*")

% Method 1: Leiame kuupsplaini punktide põhjal maksimumpunkti
% Tuletis kuupsplainist S^(3,2)(x) ehk dS/dx
dS = zeros(1, length(S3_2));
dS(1) = (S3_2(2) - S3_2(1))/hh;                                         % Ettesammuga
for i = 2:(length(S3_2)-1)
    dS(i) = (S3_2(i+1) - S3_2(i-1))/(2*hh);                             % Kesksammuga
end
dS(length(S3_2)) = ( S3_2(length(S3_2)) - S3_2(length(S3_2)-1) )/hh;    % Tahasammuga
% Leiame ekstreemumpunkti ehk S'(x) = 0 ehk see väärtus, mis 
% on nullile kõige lähemal.
% Selleks on vaja võtta igast väärtusest absoluutväärtus ja leida 
% kõige väiksem väärtus.
% Leiame väikseima väärtuse ja talle vastava indeksi:
[min_val, min_index] = min(abs(dS))
fprintf("Maksimumpunkti koordinaadid: [%.4f, %.4f]\n", xx(min_index), S3_2(min_index))

% Kontroll, kas tegemist max puntiga ehk kas S''(x) < 0
d2S = (S3_2(min_index+1) - 2*S3_2(min_index) + S3_2(min_index-1))/(h*h)
if (d2S < 0.0) 
    fprintf("Tegemist on maksimumpunktiga.\n")
else
    fprintf("Tegemist on miinimumpunktiga.\n") 
end

% Method 2: maksimumpunkti leidmine interp1()-ga
% Lähendame tuletise funktsiooni lineaarsplainiga S^(1,0)(x)
S1_0 = interp1(x, dy, x, "linear")
hold on
plot(x, S1_0, "--")
yline(0)    % punkt kus S^(1_0) läbib x-telge on funktsiooni max punkt

% Tuletis lineaarsplainist
tul = @(t) interp1(x, dy, t, "linear");
x_stats = fsolve(tul, -0.7)

[y(1) y(length(x)) interp1(x, y, x_stats, "spline")]
% max punkt
[x_stats interp1(x, y, x_stats, "spline")]


fprintf("\nÜlesanne 4.\n");

t = [0.0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6];
x = [3.1 3.8 4.9 6.5 8.4 9.9 11 11.5 11.7 11.8];
h = t(2) - t(1);

N = 3; % kuupolünoom
Vtemp = vander(t);
k = ones(1, length(t));
V = Vtemp(:, (length(Vtemp) - N):length(Vtemp)) .* k';   %V - weigthed Vandermonde matrix; N - desired polynomial degree; k - weights
C = linsolve(V, x')

p3 = @(T) C(1).*T.^3 + C(2).*T.^2 + C(3).*T + C(4);

figure(4)
plot(t, x, "o")
grid on
hold on
fplot(p3, [min(t), max(t)])

syms T
d2f = matlabFunction(diff(p3, T, 2))
hold on
fplot(d2f, [min(t), max(t)])
yline(0)
t0 = fsolve(d2f, 1.4)
x0 = p3(t0)
