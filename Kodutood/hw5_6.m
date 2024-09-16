clc
clear

fprintf("Ülesanne 1.\n")
disp("||x||_1 = |x_1| + |x_2| + ... + |x_m|")
disp("||x||_2 = sqrt(x_1^2 + x_2^2 + ... + x_m^2)")
disp("||x||_inf = max{|x_1|; |x_2|; ...; |x_m|}")

fprintf("\nÜlesanne 2.")

A = [ 8.0  7.0 -2.0
     -6.0  4.0  3.0
      4.0 -1.0  5.0];

norm(A, 1)
norm(A, inf)

fprintf("\nÜlesanne 3.")
F = @(z) [ z(1)^3 + z(2) - 1, z(2)^3 - z(1) + 1];
x0 = [0, 0];
z_lahend = fsolve(F, x0);
disp(z_lahend)

fprintf("\nÜlesanne 4.")
figure(1)
fimplicit(@(x,y) y + 5.*x.*y - x.^2)
hold on
fimplicit(@(x,y) y + x.^2 - x - 0.75)
legend("y + 5*x*y - x^2", "y + x^2 - x - 0.75")
grid on
syms x y
% Long way
%{
df1dx = diff(y + 5*x*y - x^2, x)
df1dy = diff(y + 5*x*y - x^2, y)
df2dx = diff(y + x^2 - x - 0.75, x)
df2dy = diff(y + x^2 - x - 0.75, y)
J = [ df1dx  df1dy
      df2dx  df2dy]
%}

% Short way
J = jacobian([y + 5*x*y - x^2, y + x^2 - x - 0.75],[x, y]);
disp(J)

F = @(z) [z(2) + 5*z(1)*z(2) - z(1)^2  ;
          z(2) + z(1)^2 - z(1) - 0.75 ];
dF = @(z) [5*z(2) - 2*z(1), 5*z(1) + 1 ;
           2*z(1) - 1     ,          1];


initial_values = {[-0.2; 0.5], [1.3; 0.3], [-0.5; -0.1]};
epsilon = 1e-6;
max_iter = 100;

disp("Newtoni meetod:")
[~, ~] = newton_method(F, dF, initial_values, epsilon, max_iter);
disp("Modifitseeritud Newtoni meetod:")
[~, ~] = modified_newton(F, dF, initial_values, epsilon, max_iter);

% NB! Reminder!
%x = g1(x,y)
%y = g2(x,y)
%disp("Harilik iteratsioonimeetod")
%g1 = @(x, y) (x^2 - y)/(5*y);      % Leiab ristumispunkti kui xy_vec = [-0.2; 0.5]
%g2 = @(x, y) 0.75 - x*x + x;       % Leiab ristumispunkti kui xy_vec = [-0.2; 0.5]

% Needs to be flipped!!!
%g2 = @(x,y) x^2/(1+5*x);           % Tuletised sobisid kohale [1.3; 0.3]
%g1 = @(x,y) sqrt(0.72 + x - y);    % Tuletised sobisid kohale [1.3; 0.3]

%g2 = @(x,y) x^2/(1+5*x);           % Tuletised sobivad kohale [-0.5, -0.1]
%g1 = @(x,y) (0.75 - y)/(x - 1);    % Tuletised sobivad kohale [-0.5, -0.1]
%g1 = @(x,y) nthroot(0.75 - y + x, 2);  % Tuletis sobis samale kohale, kuid koondub punkti [1.372065, 0.239502]


disp("Harilik iteratsioonimeetod:");
% Definitsioon funktsioonipaaride jaoks
g1_funcs = {@(x, y) (x^2 - y)/(5*y), @(x, y) sqrt(0.72 + x - y), @(x, y) (0.75 - y)/(x - 1)};
g2_funcs = {@(x, y) 0.75 - x*x + x,  @(x, y) x^2/(1+5*x),        @(x, y) x^2/(1+5*x)};

% Algväärtused
max_iter = 1000;
epsilon = 1E-6;
xy_vec = {[-0.2; 0.5], [1.3; 0.3], [-0.5; -0.1]};

% Iteratsioon läbi kõik funktsioonipaarid
for i = 1:length(g1_funcs)
    g1 = g1_funcs{i};
    g2 = g2_funcs{i};
    xy = xy_vec{i};
    
    fprintf('Set %d: ', i);
    [~, ~, ~] = him(g1, g2, xy, epsilon, max_iter);
end

fprintf("\nÜlesanne 5.")
figure(2)
fimplicit(@(x,y) x.*y-2)
hold on
fimplicit(@(x,y) -x+3.*y+1)
legend("x*y - 2", "-x + 3*y + 1")
grid on
initial_values = {[-2; 1], [3; 0.6]};
syms x y
F = @(z) [z(1)*z(2) - 2.0; -z(1) + 3*z(2) + 1.0];
jacobian([x*y-2, -x+3*y+1], [x,y])
dF = @(z) [ z(2), z(1); -1.0, 3.0];

epsilon = 1e-6;
max_iter = 100;

disp("Newtoni meetod:")
[~, ~] = newton_method(F, dF, initial_values, epsilon, max_iter);
disp("Modifitseeritud Newtoni meetod:")
[~, ~] = modified_newton(F, dF, initial_values, epsilon, max_iter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z, count] = newton_method(F, dF, initial_values, epsilon, max_iter)
    % F - funktsioonide veeruvektor, mille nullkohta otsime
    % dF - jakobiaani maatriks
    % initial_values - alglähendide paarid / ristumispunktid
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

function [z, count] = modified_newton(F, dF, initial_values, epsilon, max_iter)
    % F - funktsioonide veeruvektor, mille nullkohta otsime
    % dF - jakobiaani maatriks
    % initial_values - alglähendide paarid / ristumispunktid
    % epsilon - täpsus
    % max_iter - maksimaalne iteratsioonide arv

    % Käime läbi iga arvupaari
    for i = 1:length(initial_values)
        z = initial_values{i};
        z0 = z;
        count = 0;
          
        while norm(F(z), 1) > epsilon && count < max_iter
            z = z - dF(z0)^(-1) * F(z); % modifitseeritud Newtoni meetodi eeskiri
            count = count + 1;
        end
          
        fprintf('Set %d: Iterations = %d, Solution = [%.6f, %.6f]\n', i, count, z(1), z(2));
    end
end

% J. Janno lk. 70
function [x, y, count] = him(g1, g2, xy_vec, epsilon, max_iter)
    x = xy_vec(1);
    y = xy_vec(2);
    xvana = xy_vec(1) - 1;
    yvana = xy_vec(2) - 1;
    count = 0;
    while max(abs(x-xvana), abs(y-yvana)) >= epsilon && count < max_iter
        xvana = x;
        yvana = y;
        x = g1(x,y);
        y = g2(xvana,y);
        count = count + 1;
    end
    if count == max_iter
        warning('Maximum number of iterations reached without convergence.');
    end
    fprintf('Iterations = %d, Solution = [%.6f, %.6f]\n', count, x, y);
end
