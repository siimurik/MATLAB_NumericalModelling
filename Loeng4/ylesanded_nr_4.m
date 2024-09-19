clc
clear

fprintf("\nÜlesanne 1.")

figure(1)
fplot(@(x)x.^6 - x - 1 )
xlim([-1.5, 1.5])
ylim([-4, 4])
xline(0)
yline(0)
grid on

syms x
%g = @(x) x.^6 - 1.0;                   % ei sobi, hajub
g1 = @(x) -nthroot(x + 1.0, 6);   % plus ja miinus peab ka ees olema
g2 = @(x) +nthroot(x + 1.0, 6);   % plus ja miinus peab ka ees olema
%dg = matlabFunction(diff(x.^6 - 1.0)); % ei sobi, hajub
dg = matlabFunction(diff(nthroot(x + 1.0, 6))); % vahet pole kas - või + panna. sama tulemus
x0 = [-0.75, 1.2];
disp(abs(dg(x0)))
disp(abs(dg(x0)) < 1)

%{
LISA:
syms x
tuletis = diff(g, x)

syms tul(x)
tul(x) = tuletis

% kui ei meeldi harilik murd siis vpa() teisendab kümnendmurruks 32 komakohaga
vpa(abs(tul(-0.75)))    
vpa(abs(tul(1.2)))    
%}
fprintf("Harilik iteratsioonimeetod:\n")
epsilon = 1.0E-6;
max_iter = 1000;
[~,~] = him(g1, x0(1), epsilon, max_iter);
[~,~] = him(g2, x0(2), epsilon, max_iter);

fprintf("Newtoni meetod:\n")
f = @(x) x.^6 - x - 1.0;
df = @(x) 6.*x.^5 - 1.0;

[~,~] = newton_method(f, df, x0(1), epsilon, max_iter);
[~,~] = newton_method(f, df, x0(2), epsilon, max_iter);

fprintf("\nÜlesanne 3.")
fprintf("Harilik iteratsioonimeetod:\n")
[~,~] = him(g2, 0.0, epsilon, max_iter); % g1-korral hajub
fprintf("Newtoni meetod:\n")
[~,~] = newton_method(f, df, 0.0, epsilon, max_iter);

fprintf("\nÜlesanne 2.")
clear
figure(2)
fplot(@(x) x.^3 - 6.*x.^2 + 3.*atan(x) + 3)
xline(0); yline(0)
xlim([-2, 7])
ylim([-30, 30])
grid on

fprintf("Harilik iteratsioonimeetod:\n")
syms x
%g = @(x) tan((-x.^3 + 6.*x.^2 - 3)./3.0);
g1 = @(x) nthroot(6.*x.^2 - 3.*atan(x) - 3.0, 3);
dg1 = matlabFunction(diff( nthroot(6.*x.^2 - 3.*atan(x) - 3.0, 3) ));
g2 = @(x) sqrt((x.^3 + 3.*atan(x) + 3)./6.0);
dg2 = matlabFunction(diff( sqrt((x.^3 + 3.*atan(x) + 3)/6.0) ));
x0 = [-0.5, 1.0, 5.8];
abs(dg1(x0))    % sobib 5.8 jaoks
abs(dg2(x0))    % sobib -0.5, 1.0 jaoks

epsilon = 1.0E-6;
max_iter = 1000;
g2_neg = @(x) -sqrt((x.^3 + 3.*atan(x) + 3)/6.0);
[~, ~] = him(g2_neg, x0(1), epsilon, max_iter);
[~, ~] = him(g2,     x0(2), epsilon, max_iter);
[~, ~] = him(g1,     x0(3), epsilon, max_iter);

fprintf("Newtoni meetod:\n")
f = @(x) x.^3 - 6.*x^2 + 3.*atan(x) + 3;
df = matlabFunction(diff( x.^3 - 6.*x^2 + 3.*atan(x) + 3 ));

for i = 1:length(x0)
    [~, ~] = newton_method(f, df, x0(i), epsilon, max_iter);
end

fprintf("Modifitseeritud Newtoni meetod:\n")
for i = 1:length(x0)
    [~, ~] = modified_newton(f, df, x0(i), epsilon, max_iter);
end

fprintf("\nÜlesanne 4.")
f = @(x) x.^5 - x.^2 - 2.0;
figure(3)
fplot(f)
grid on
xlim([-2, 2.5])
ylim([-30, 30])
xline(0)
yline(0)

[~, ~] = secant_method(f, 1.0, 1.5, epsilon, max_iter);


function [x_uus, count] = him(g, x_algne, tol, max_iter)
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
    while abs(x_uus - x_vana) >= tol && count < max_iter
        x_vana = x_uus;
        x_uus  = g(x_vana); % h.i.m eeskiri
        count  = count + 1;
    end

    % Kontroll, kas maksimaalne iteratsioonide arv ületati
    if count >= max_iter
        warning('Maksimaalne iteratsioonide arv ületatud.');
    end

    fprintf('Iterations = %d, Solution = [%.6f, %.6f]\n', count, x_uus, 0.0);
end

function [x, count] = newton_method(f, df, x0, epsilon, max_iter)
    % f - funktsioon, mille nullkohta otsime
    % df - funktsiooni tuletis
    % x0 - alglähend
    % epsilon - täpsus
    % max_iter - maksimaalne iteratsioonide arv
    x1 = x0;
    count = 0;
    while abs(f(x1)) > epsilon && count < max_iter
        x0 = x1; % Uuenda x0 väärtusega x1
        x1 = x0 - f(x0)/df(x0); % Newtoni meetodi
        count = count + 1;
    end
    x = x1;

    % Kontroll, kas maksimaalne iteratsioonide arv ületati
    if count >= max_iter
        warning('Maksimaalne iteratsioonide arv ületatud.');
    end

    fprintf('Iterations = %d, Solution = [%.6f, %.6f]\n', count, x, 0.0);
end

function [x, count] = modified_newton(f, df, x0, eps, max_iter)
    % f - funktsioon, mille nullkohta otsime
    % df - funktsiooni tuletis
    % x0 - alglähend
    % eps - täpsus
    % max_iter - maksimaalne iteratsioonide arv

    % Algväärtustamine
    x = x0;
    count = 0;
    df_x0 = df(x0); % Tuletise väärtus alglähendis

    % Iteratiivne protsess
    while abs(f(x)) > eps && count < max_iter
        x0 = x;
        x  = x0 - f(x0) / df_x0;
        count = count + 1;
    end

    % Kui maksimaalne iteratsioonide arv on saavutatud
    if count >= max_iter
        warning('Maksimaalne iteratsioonide arv ületatud.');
    end

    fprintf('Iterations = %d, Solution = [%.6f, %.6f]\n', count, x, 0.0);
end

function [x, count] = secant_method(f, x0, x1, eps, max_iter)
    % f - funktsioon, mille nullkohta otsime
    % x0 - esimene alglähend
    % x1 - teine alglähend
    % eps - täpsus
    % max_iter - maksimaalne iteratsioonide arv

    count = 0;
    while (abs(f(x1))) > eps && count < max_iter
        x2 = x1 - f(x1)*(x1 - x0)/(f(x1) - f(x0));
        x0 = x1;
        x1 = x2;
        count = count + 1;
    end

    % Kontroll, kas maksimaalne iteratsioonide arv ületati
    if count >= max_iter
        warning('Maksimaalne iteratsioonide arv ületatud.');
    end
    x = x1;

    fprintf('Iterations = %d, Solution = [%.6f, %.6f]\n', count, x, 0.0);
end