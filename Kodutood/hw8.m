clc 
clear

fprintf("Ülesanne 1.\n")
fprintf("a) Funktsiooni graafiline kujutis\n")
figure(1)
fplot(@(x) x.^5 - 2.*x.^3 + 2.*x - 2)
xlim([-2, 2.5])
ylim([-20, 20])
grid on

g = @(x) nthroot(2.*x.^3 - 2.*x + 2, 5);
x_init  = 1.1;
tol      = 1.0E-6;
max_iter = 1000;
[~, ~] = him(g, x_init, tol, max_iter);

fprintf("\nÜlesanne 2.")
c = [-4, 3, 6, 1, 0, -7];
fprintf("||c||_2   = %.6f\n", norm(c, 2))
fprintf("||c||_1   = %.6f\n", norm(c, 1))
fprintf("||c||_inf = %.6f\n", norm(c, inf))

fprintf("\nÜlesanne 3.")

A = [-5 3; 2 -1];
B = [4 -2; 7  0];

C = A .* B;
disp(C)
fprintf("B(2,1) = %.4f\n", B(2, 1))

fprintf("\nÜlesanne 4.")
y = @(x) x.^3 - 2.*x.^2 + 7;
fprintf("y(5) = %f\n", y(5))
x = [0, 1, 4, 3];
y(x)

fprintf("\nÜlesanne 5.")
F = @(z) [  z(1)^2 - 2*z(2) + 4;
            z(2)^2 - z(1)^4 + 4*z(1) - 1 ];
syms x y
jacobian([x^2 - 2*y + 4, y^2 - x^4 + 4*x - 1] , [x,y])
dF = @(z) [      2*z(1),     -2;
           4 - 4*z(1)^3, 2*z(2)];

initial_values = {[-1, 2.5]};
epsilon  = 1.0E-6;
max_iter = 1000;
[~, ~] = newton_method(F, dF, initial_values, epsilon, max_iter);

fprintf("\nÜlesanne 6.")
f = @(x,y) cos(x.*y) + 3.*exp(sin(7.*x)) - nthroot(2.*x.^3.*y.^4 + 4, 5);
figure(2)
fimplicit(f)
grid on

figure(3)
x = linspace(-3, 3);
y = linspace(-3, 3);
[X, Y] = meshgrid(x,y);
Z = f(X, Y);
surf(X, Y, Z)

syms x y
F = cos(x*y) + 3*exp(sin(7*x)) - nthroot(2*x^3*y^4 + 4, 5);
dFdx = diff(F, x)
dFdy = diff(F, y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_uus, count] = him(g, x_init, tol, max_iter)
    % Leiab lahendi hariliku iteratsioonimeetodiga
    % g - iteratsioonifunktsioon
    % x_init - algväärtus
    % tol - tolerants
    % max_iter - maksimaalne iteratsioonide arv

    % Algväärtused
    x_uus  = x_init;
    x_vana = x_init - 1; % Algväärtus, mis erineb x_uus-st
    count  = 0;

    % Iteratsiooniprotsess
    while abs(x_uus - x_vana) >= tol && count < max_iter
        x_vana = x_uus;
        x_uus  = g(x_vana);
        count  = count + 1;
    end

    % Kontroll, kas maksimaalne iteratsioonide arv ületati
    if count >= max_iter
        warning('Maksimaalne iteratsioonide arv ületatud.');
    end

    fprintf('Iterations = %d, Solution = [%.6f, %.6f]\n', count, x_uus, 0.0);
end

function [z, count] = newton_method(F, dF, initial_values, epsilon, max_iter)
    % F - funktsioonide veeruvektor, mille nullkohta otsime
    % dF - jakobiaani maatriks
    % initial_values - alglähendide paarid / ristumispunktid. initial_values={[...],...,[...]}
    % epsilon - täpsus
    % max_iter - maksimaalne iteratsioonide arv

    % Käime läbi iga arvupaari
    for i = 1:length(initial_values)
        z = initial_values{i};
        count = 0;
          
        while norm(F(z), 1) > epsilon && count < max_iter
            z0 = z;
            z  = z0 - dF(z0)^(-1) * F(z0); % Newtoni meetodi eeskiri
            count = count + 1;
        end
          
        fprintf('Set %d: Iterations = %d, Solution = [%.6f, %.6f]\n', i, count, z(1), z(2));
    end
end