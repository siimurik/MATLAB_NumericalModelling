clc
clear

fprintf("Ülesanne 1.\n")
fprintf("a) p(x) väärtus kui x=0 ja x=6\n")
p = @(x) 4.*x.^5 - 2.*x.^4 + x.^3 + 7.*x - 9.0;
fprintf("p(x = %d) = %.4f\n", 0, p(0))
fprintf("p(x = %d) = %.4f\n", 6, p(6))

r = [4 -2 1 0 7 -9];
% p(0)
polyval(r, 0)
polyval(r, 6)

% polünoom p(x) graafikul -1 < x < 8
fprintf("b) p(x) graafik piirkonnas -1 < x < 8\n")
% 1. variant
p = @(x) 4.*x.^5 - 2.*x.^4 + x.^3 + 7.*x - 9.0;
figure(1)
fplot(p, [-1, 8])
grid on
hold off

% 2. variant
x = -1:0.01:8;
y = polyval(r, x);
figure(2)
plot(x, y)  % võib vajadusel lisada ka "-o"

fprintf("\nÜlesanne 2.\n")
clear
p = [1.0 -1.0 -12.0];
roots(p)    % polünoomi juured
% Polynomial differentation - polyder()

fprintf("\nÜlesanne 3.\n")
juured = [-1.0, 3.0, 2.0, 5.0];
p = poly(juured)    % polünoomikordajad
x = -2:0.01:7;      % vahemik suurem kui etteantud juured
y = polyval(p, x);
figure(3)
plot(x,y)
grid on
hold on
x0LinePoints = zeros(1, length(juured));    % juured graafiliselt
plot(juured, x0LinePoints, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
hold off

fprintf("\nÜlesanne 4.\n")
p1 = [1.0 -2.0, 1.0];
p2 = [1.0  2.0];
conv(p1, p2)    % polünoomide korrutamine
fprintf("Konvolutsiooni käsk annab meile %d. astme polünoomi.\n", length(conv(p1, p2)))

fprintf("\nÜlesanne 5.\n")
x = [-1 0.1 0.5 3 4 6.3  7  9 14 21 ];
y = [-2 0.4 0.7 2 4 3.6 3.8 6 -1 12 ];
figure(4)
plot(x, y, 'o')
grid on
hold on
% lähendame erinevate polünoomidega
% polünoomid joonisele illustratiivsuse pärast
% POLYFIT Fit polynomial to data.
%   P = POLYFIT(X,Y,N) finds the coefficients of a polynomial P(X) of
%   degree N that fits the data Y best in a least-squares sense.
% 1. astme polünoom
p1 = polyfit(x, y, 1)   % ehk polünoomi kordajad
% p1(x) = 0.4172*x + 0.2839
xx = [-3:0.01:23];
y1 = polyval(p1, xx);   % leidsime p1(xx)
plot(xx, y1)
% 2. astme polünoom
p2 = polyfit(x, y, 2)   % ehk polünoomi kordajad
y2 = polyval(p2, xx);   % leidsime p2(xx)
plot(xx, y2)
% 3. astme polünoom
p3 = polyfit(x, y, 3)   % ehk polünoomi kordajad
y3 = polyval(p3, xx);   % leidsime p3(xx)
plot(xx, y3)
% 6. astme polünoom
p6 = polyfit(x, y, 6)   % ehk polünoomi kordajad
y6 = polyval(p6, xx);   % leidsime p6(xx)
plot(xx, y6)
hold off

% polünoomi väärtus kui x = 3.5
x_val = 3.5;
y_val_p2 = polyval(p2, x_val);
y_val_p6 = polyval(p6, x_val);
disp(['Vastus 2. astme polünoomi korral kui x = 3.5: ', num2str(y_val_p2)]);
disp(['Vastus 6. astme polünoomi korral kui x = 3.5: ', num2str(y_val_p6)]);

fprintf("\nÜlesanne 7.\n")
x_aasta = [2014 2015 2016 2017 2018 2019 2020 2021];
y_oC    = [22.2 19.5 18.0 17.9 25.2 18.6 18.7 29.3];
figure(5)
plot(x_aasta, y_oC, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
grid on
hold on
fprintf("Lähendamine 1. astme polünoomiga:\n")
p1 = polyfit(x_aasta, y_oC, 1)   % ehk polünoomi kordajad
xdata = linspace(min(x_aasta)-1, max(x_aasta)+1, 1000);
y1 = polyval(p1, xdata);
plot(xdata, y1)
fprintf("Lähendamine 4. astme polünoomiga:\n")
p4 = polyfit(x_aasta, y_oC, 4)   % ehk polünoomi kordajad
y4 = polyval(p4, xdata);
plot(xdata, y4)

x_prognoos = 2024;
y_prog_p1 = polyval(p1, x_prognoos);
y_prog_p4 = polyval(p4, x_prognoos);
fprintf("Prognoos 1. astme polünoomi kohaselt on 2024. aastal temp. %.4f C.\n", y_prog_p1)
fprintf("Prognoos 4. astme polünoomi kohaselt on 2024. aastal temp. %.4f C.\n", y_prog_p4)

fprintf("\nÜlesanne 6.\n")
x = [-1.0  0.0 1.0  2.0];
y = [ 3.0 -4.0 5.0 -6.0];
figure(6)
plot(x, y, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')

syms X
L_3_0 = expand(((X-0)*(X-1)*(X-2))/((-1-0)*(-1-1)*(-1-2)))
L_3_1 = expand(((X-(-1))*(X-1)*(X-2))/((0-(-1))*(0-1)*(0-2)))
L_3_2 = expand(((X-(-1))*(X-0)*(X-2))/((1-(-1))*(1-0)*(1-2)))
L_3_3 = expand(((X-(-1))*(X-0)*(X-1))/((2-(-1))*(2-0)*(2-1)))

lagrp = matlabFunction(L_3_0*y(1) + L_3_1*y(2) + L_3_2*y(3) + L_3_3*y(4))

% interpol tingimused
lagrp(x)    % samad tulemused mis y-vektoris
hold on 
fplot(lagrp)

lagrp(0.5)
