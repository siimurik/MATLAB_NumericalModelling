clc
clear
format short
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
title('Ül. 2 Funkstiooni graafik');
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
title("Ül.3 Funktsiooni graafik")
yline(0, 'k-'); % x-telje joon
xlim([-0.5, 5.])
ylim([-10, 10])
grid on

fprintf("\nb) fsolve() lahend\n")
x0 = [0, 2, 3.5]; % Punktid, mis visuaalselt umbes läbivad x-telge
x  = fsolve(y,x0);
disp(x)
% Lisan nullkohad eelmisele graafikule
hold on
plot(x, [0,0,0], 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
hold off

fprintf("\nc) Lahend hariliku iteratsioonmeetodiga")
%   2.*x.^3 - 11.7.*x.^2 + 17.7.*x - 5 = 0
%g1 = @(x) 5.0./(2.*x.^2 - 11.7.*x + 17.7);                 
g1 = @(x) (5.0 - 2.*x.^3 + 11.7.*x.^2)./17.7;               % Leiab esimese nullkoha. 
g2 = @(x) sqrt((2.*x.^3 + 17.7.*x - 5.0)./11.7);            % Leiab teise nullkoha.
g3 = @(x) nthroot((5.0 + 11.7.*x.^2 - 17.7.*x)./2.0, 3);    % Leiab kolmanda nullkoha.
G = {g1, g2, g3};
x0 = [0, 2, 3.5]; % Punktid, mis visuaalselt umbes läbivad x-telge
%x_algne = 0.0;  
tol = 1.0E-6;
max_iter = 1000000;
%[x_uus, count] = harilik_iter(g3, x_algne, tol, max_iter);
for i = 1:length(x0)
    [~, ~] = him(y, G{i}, x0(i), tol, max_iter);
end


fprintf("\nÜlesanne 4.")
alpha = 4.0;
beta  = 5.0;
gamma = 1.0;

y = @(x) cos(x) + (-1)^gamma.*0.1.*alpha.*x - 0.001*beta;
figure(3)
fplot(y)
title("Ül.4 Funktsiooni graafik")
yline(0, 'k-'); % x-telje joon
grid on
hold off

x0 = 1.0; % Punkt, mis visuaalselt umbes läbib x-telge
x  = fsolve(y,x0);
disp(x)

fprintf("\nÜlesanne 5.")
a = 4.0;
b = 5.0;
c = 1.0;
figure(4)
y = @(x) a./10.*x - b./40 - (-1).^c.*sin(x);
yline(0, 'k-'); % x-telje joon
fplot(y)
title("Ül.5 Esiaglne funktsioon")
grid on
hold off

% Defineerides 
%   a/10*x - b/40 - (-1)^c*sin(x) = 0
%   a/10*x = b/40 + (-1)^c*sin(x)
%   x = 10/a*(b/40 + (-1)^c*sin(x))
%g = @(x) 10./a*(b/40. + (-1.0)^c*sin(x));

figure(5)
fplot(@(x) asin(1./(-1).^c) .* ((a./10.).*x - b./40.))      
%fplot(@(x) 10./a.*(b./40. + (-1.0).^c*sin(x)), [-10, 10])  
grid on                                                     
yline(0, 'k-'); % x-telje joon
hold off

g = @(x) asin(1./((-1).^c) .* ((a./10.).*x - b./40.));  % Leiab lahendi
%g = @(x) 10./a*(b/40. + (-1.0)^c*sin(x));              % Ei leia lahendit; tingimus |g'(x)| < 1 pole rahuldatud
x_algne = 0.2;
tol = 1.0E-15;
max_iter = 100000;
[x_uus, count] = harilik_iter(g, x_algne, tol, max_iter);
fprintf("x_uus = %.4f\n", x_uus)
fprintf("count = %d\n", count)
%[x0, k] = him(y, g, x_algne, tol, max_iter);
%fprintf("x_0 = %.15f\n", x0)
%fprintf("k = %d\n", k)
% Lisan uue leitud nullkoha originaalse graafiku joonisele
figure(4)
hold on
plot(x_uus, 0.0, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
hold off

% Vastuse kontrollimine
fprintf("\nKontroll fsolve() funkstiooni kasutades:")
x0 = 0.2;
xf = fsolve(@(x) a./10.*x - b./40 - (-1).^c.*sin(x), x0);
disp(xf)

function [x_uus, count] = harilik_iter(g, x_algne, tol, max_iter)
    % Leiab lahendi hariliku iteratsioonimeetodiga
    % g - iteratsioonifunktsioon
    % x_algne - algväärtus
    % tol - tolerants
    % max_iter - maksimaalne iteratsioonide arv

    % Algväärtused
    x_uus  = x_algne;
    x_vana = x_algne - 1; % Algväärtus, mis erineb x_uus-st
    count  = 0;

    % Iteratsiooniprotsess
    while norm(abs(x_uus - x_vana), 2) >= tol && count < max_iter
        x_vana = x_uus;
        x_uus  = g(x_vana);
        count  = count + 1;
    end

    % Kontroll, kas maksimaalne iteratsioonide arv ületati
    if count >= max_iter
        warning('Maksimaalne iteratsioonide arv ületatud.');
    end
end

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

    fprintf('Iterations = %d, Solution = [%.6f]\n', count, x_uus);

    % Kontroll, kas maksimaalne iteratsioonide arv ületati
    if count >= max_iter
        warning('Maksimaalne iteratsioonide arv ületatud.');
    end
end

%{
    "Juureotsingu" funktsioon 
---------------------------------------------------
Kasutades ROOTS käsku leiab see käsk polünoomi juured, 
mille kordajad on antud vektoris c. Kui vektoril c 
on (N+1) elementi, siis on polünoom kujul
    C(1)*X^N + ... + C(N)*X + C(N+1) .
Lahendamiseks koostadakse vektorist c kaasmaatriks. 
ROOTS leiab antud kaasmaatriksi omaväärtused, mis
on antud polünoomi nullkohad.
Vt. 'kaasmaatriks' -> [Link]
https://en.wikipedia.org/wiki/Companion_matrix
---------------------------------------------------
LISA: selgitus, mida ROOTS funktsioon teeb:
Olgu antud polünoom
p(x) = 2*x^5 - 5*x^4 + 8*x^3 - 9*x^2 + 10*x - 12 .
Lisame koefitsendid c vektorisse
    >> c = [2 -5 8 -9 10 -12]
Leiame kaasmaatriksi.
    >> A = compan(c)
A =
    2.5   -4.0    4.5   -5.0    6.0
    1.0      0      0      0      0
      0    1.0      0      0      0
      0      0    1.0      0      0
      0      0      0    1.0      0

Selle käsklusega normaliseeritakse algne funktsioon
maatriksi A esimeses reas uule kujule, kus iga 
kordaja jagatakse läbi c(1) = 2 liikmega:
    p(x) = 2/2*x^5 - 5/2*x^4 + 8/2*x^3 - ...
            - 9/2*x^2 + 10/2*x - 12/2 .
Seega on funksiooni p(x) uus kuju
    p(x) = x^5 - 2.5*x^4 + 4*x^3 -...
            - 4.5*x^2 + 5*x - 6 .
Viimase sammuna on on iga kordaja vaja korrutada läbi 
(-1)-ga. Seega kui meil on praegu funktsioon kujul
    p(x) = x^5 + a_4*x^4 + a_3*x^3 +...
            +a_2*x^2 + a_1*x + a_0
Siis maatrikisse lisamisel on vormistus järgmine:
A = [-a_4 -a_3 -a_2 -a_1 -a_0
      1.0  0.0  0.0  0.0  0.0
      0.0  1.0  0.0  0.0  0.0
      0.0  0.0  1.0  0.0  0.0
      0.0  0.0  0.0  1.0  0.0]
Nagu näha, esimeses reas on polünoomi p(x) kordajad 
ning sellele järgnevatel ridadel lisatakse peadia-
gonaalist ühe sammu võrra alla ühed kuni viimase reani,
kus viimase rea ja veeru element peab olema võrdne
nulliga nagu ülejäänud maatriksi elemendid.

Vastuse esitamiseks leiame maatriksi A omaväärtused:
    >> eig(A)
ans =
   1.5182 + 0.0000i
   0.9211 + 1.2989i
   0.9211 - 1.2989i
  -0.4301 + 1.1721i
  -0.4301 - 1.1721i
%}
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
            fprintf("[%f, %f]\n", real(r(i)), 0.0)
            x_zeros = [x_zeros, real(r(i))];
        end
    end
end