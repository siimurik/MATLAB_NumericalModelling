clc
clear

% Iseseisev töö nr. 4
fprintf("Ülesanne 1.\n")
%syms x
%diff(x^4 + 2*x^2 - x - 3)
f  = @(x) x.^4 + 2.*x.^2 - x - 3;
df = @(x) 4.*x.^3 + 4.*x - 1;
figure(1)
fplot(f)
%plot_darkmode
yline(0, "k-", "LineWidth", 1.)
xlim([-3, 3])
ylim([-8, 6])
grid on

x0  = [-1.0, 1.0];
eps = 1E-3;
max_iter = 100000;
fprintf("\nTäpsusaste iga meetodi jaoks: eps = %.1e.", eps)

% Newtoni meetod
fprintf("\nLahendus Newtoni meetodiga:")
[x, count] = newton_method(f, df, x0(1), eps, max_iter);
fprintf("x_(n+1) = %.4f\n", x)
fprintf("count = %d\n",  count)
hold on; plot(x, 0, "ro", "MarkerSize", 6, "MarkerFaceColor", "r")
[x, count] = newton_method(f, df, x0(2), eps, max_iter);
fprintf("x_(n+1) = %.4f\n", x)
fprintf("count = %d\n",  count)
hold on; plot(x, 0, "ro", "MarkerSize", 6, "MarkerFaceColor", "r")
hold off

% Modifitseeritud Newtoni meetod
fprintf("\nLahendus modifitseeritud Newtoni meetodiga:")
[x, count] = modified_newton(f, df, x0(1), eps, max_iter);
fprintf("x_(n+1) = %.4f\n", x)
fprintf("count = %d\n",  count)
[x, count] = modified_newton(f, df, x0(2), eps, max_iter);
fprintf("x_(n+1) = %.4f\n", x)
fprintf("count = %d\n",  count)

% Lõikajate meetod
fprintf("\nLahendus lõikajate meetodiga:")
x0 = -1.0;  x1 = -0.5;
[x, count] = secant_method(f, x0, x1, eps, max_iter);
fprintf("x_(n+1) = %.4f\n", x)
fprintf("count = %d\n",  count)
x0 =  1.0;  x1 =  1.5;
[x, count] = secant_method(f, x0, x1, eps, max_iter);
fprintf("x_(n+1) = %.4f\n", x)
fprintf("count = %d\n",  count)

fprintf("\nÜlesanne 2.\n")
a = 4.0;
b = 5.0;
c = 1.0;
f = @(x)  a./10.*x - b./40 - (-1).^c.*sin(x);
figure(2)
fplot(f);
yline(0, "k-", "LineWidth", 1.)
grid on

syms x
deriv = diff( a/10*x - b/40 - (-1)^c*sin(x));
fprintf("Funktsiooni f tuletis on %s.\n", deriv)
df = @(x) cos(x) + 2./5;

% Algväärtused
x0  = 0.1;
eps = 1E-6;
max_iter = 100000;
fprintf("\nTäpsusaste iga meetodi jaoks: eps = %.1e.", eps)

% Newtoni meetod
fprintf("\nLahendus Newtoni meetodiga:")
[x, count] = newton_method(f, df, x0, eps, max_iter);
fprintf("x_(n+1) = %.4f\n", x)
fprintf("count = %d\n",  count)

% Modifitseeritud Newtoni meetod
fprintf("\nLahendus modifitseeritud Newtoni meetodiga:")
[x, count] = modified_newton(f, df, x0, eps, max_iter);
fprintf("x_(n+1) = %.4f\n",  x)
fprintf("count = %d\n",  count)

hold on; plot(x, 0, "ro", "MarkerSize", 6, "MarkerFaceColor", "r")
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, count] = newton_method(f, df, x0, eps, max_iter)
    % f - funktsioon, mille nullkohta otsime
    % df - funktsiooni tuletis
    % x0 - alglähend
    % eps - täpsus
    % max_iter - maksimaalne iteratsioonide arv
    x1 = x0;
    count = 0;
    while abs(f(x1)) > eps && count < max_iter
        x0 = x1; % Uuenda x0 väärtusega x1
        x1 = x0 - f(x0)/df(x0);
        count = count + 1;
    end
    x = x1;

    % Kontroll, kas maksimaalne iteratsioonide arv ületati
    if count >= max_iter
        warning('Maksimaalne iteratsioonide arv ületatud.');
    end
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
    while (abs(f(x))) > eps && count < max_iter
        x0 = x;
        x  = x0 - f(x0) / df_x0;
        count = count + 1;
    end

    % Kui maksimaalne iteratsioonide arv on saavutatud
    if count == max_iter
        warning('Maksimaalne iteratsioonide arv saavutatud.');
    end
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
end

