clc
clear

fprintf("Ülesanne 1.\n")
x = 3: 0.1: 3.6;
y = [5.26, 5.32, 5.41, 5.37, 5.31, 5.25, 5.16];
h = 0.1;
figure(1)
plot(x,y,'o')
grid on
xlabel('x-telg')
ylabel('y-telg')
title('Andmed')
hold off


% Preallocation for speed
tuletis = zeros(1, length(y));
%-------------------------------------------
% a) tuletis esimeses sõlmes
% kasutame diferentsvalemit sammuga ette
% sammu pikkus h=0.1
tuletis(1) = (y(2)-y(1))/h;
% b) tuletis sisesõlmedes
% kasutame keskmistatud diferentsvalemit -
% tuletise leidmisel sõlmes x_i kasutame
% eelmist sõlme x_(i-1) ja järgmist sõlme x_(i+1)
for i=2:(length(y)-1)
    tuletis(i) = (y(i+1)-y(i-1))/(2*h);
end
% ac) tuletis viimases sõlmes
% kasutame diferentsvalemit sammuga taha
% sammu pikkus h=0.1
tuletis(length(y)) = (y(length(y))-y(length(y)-1))/h;
disp(tuletis)

fprintf("Ülesanne 2.\n")
x = [-2   0.0 1.5 3.2  4.9  6.1  8.1  9.2];
y = [-2.8 2.1 5.5 9.7 13.5 16.2 21.1 23.8];
figure(2)
plot(x,y,'o')
grid on
xlabel('x')
ylabel('y')
title('Andmed')

% a) lähendame funktsiooni y(x) lineaarfunktsiooniga
% vähimruutude mõttes kaalusid mitte arvestades,
% st kaalud on võrdsed
kordajad = polyfit(x, y, 1)
% y(x) ~ Phi(x) = 2.3598*x  +  1.9933
phi = @(x) kordajad(1).*x + kordajad(2);
hold on
fplot(phi,'g',[-3, 24])
% b) leiame kuupfunktsiooni 1. järku tuletise, mis annab meile tõusu
syms x;
f2 = matlabFunction(diff(phi,x,1))

fprintf("Ülesanne 3.\n")
t = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0];
s = [0.0, 0.7, 1.8, 3.4, 5.1,  6.3,  7.3,  8.0,  8.4];

fprintf("Kiiruse definitsioon: v = ds/dt\n")
dt = t(2)-t(1);
v = zeros(1, length(s));    % Et oleks kiirem
v(1) = (s(2)-s(1))/dt;                              % ettesamm
for i = 2:(length(s)-1)
    v(i) = (s(i+1) - s(i-1))/(2*dt);                % kesksamm
end
v(length(s)) = (s(length(s))-s(length(s)-1))/dt;    % tahasamm
% Leia mitmes element on t vektroris võrdne 10-ga ning 
% leia selle indeksiga element kiiruse vektorist.
for i = 1:length(t)
    if t(i) == 10.0
        fprintf("Indeks, kus t = 10 sekundit: t(%d).\n", i)
        fprintf("Seega on kiirus ajahetkel t=10 võimalik leida järgmiselt: v(%d) = %.4f m/s.\n", i, v(i))
        break
    end
end

fprintf("\nKiirenduse definitsioon: a = dv/dt")
a = zeros(1, length(v));
a(1) = (v(2) - v(1))/dt;                            % ettesamm
for i = 2:(length(v)-1)
    a(i) = (v(i+1) - v(i-1))/(2*dt);                % kesksamm
end
a(length(v)) = (v(length(v)) - v(length(v)-1))/dt;  % tahasamm
for i = 1:length(t)
    if t(i) == 10.0
        fprintf("Indeks, kus t = 10 sekundit: t(%d).\n", i)
        fprintf("Seega on kiirendus ajahetkel t=10 võimalik leida järgmiselt: a(%d) = %.4f m/s^2.\n", i, a(i))
        break
    end
end