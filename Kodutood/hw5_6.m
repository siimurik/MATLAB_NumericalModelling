clc
clear


fprintf("\nÜlesanne 2.\n")

A = [ 8.0  7.0 -2.0
     -6.0  4.0  3.0
      4.0 -1.0  5.0];

norm(A, 1)
norm(A, inf)

fprintf("\nÜlesanne 3.\n")
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
%syms x y
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
%J = jacobian([y + 5*x*y - x^2, y + x^2 - x - 0.75],[x, y]);
%disp(J)

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

% Fixed-point iteration functions for the system of equations
g1 = @(x) 0.75 - x^2 + x;          % Rearranged second equation: y in terms of x
g2 = @(x, y) x^2 / (1 + 5*x);      % Rearranged first equation: y in terms of x and y

% Initial guesses for x and y
x0 = 1.3;
y0 = g1(x0);  % Initial y guess from g1

% Tolerance and maximum iterations
tol = 1e-6;
max_iter = 1000;

% Fixed-point iteration method for solving the system of equations
[solution_x, solution_y, iterations] = fixed_point_iteration_system(g1, g2, x0, y0, tol, max_iter);
disp(['Solution: x = ', num2str(solution_x), ', y = ', num2str(solution_y)]);
disp(['Iterations: ', num2str(iterations)]);

%{
G = @(z) [sqrt(z(2)+5*z(1)*z(2)); sqrt(z(1) - z(2) + 0.75)];
jacobian([sqrt(y+5*x*y), sqrt(x - y + 0.75)],[x, y])
dG = @(z) [(5*z(2))/(2*(z(2) + 5*z(1)*z(2))^(1/2)), (5*z(1) + 1)/(2*(z(2) + 5*z(1)*z(2))^(1/2));
            1/(2*(z(1) - z(2) + 3/4)^(1/2)),      -1/(2*(z(1) - z(2) + 3/4)^(1/2))            ];
[z, count] = him(G, dG, initial_values, epsilon, max_iter);
%}

fprintf("\nÜlesanne 5.")
figure(2)
fimplicit(@(x,y) x.*y-2)
hold on
fimplicit(@(x,y) -x+3.*y+1)
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

% Fixed-point iteration method for the system
function [x, y, count] = fixed_point_iteration_system(g1, g2, x0, y0, tol, max_iter)
    x = x0;
    y = y0;
    count = 0;

    for i = 1:max_iter
        % Update y using g1 (depends only on x)
        y_new = g1(x);
        
        % Update x using g2 (depends on both x and y)
        x_new = sqrt(g2(x, y_new));  % Take the square root since g2 gives x^2

        % Check for convergence
        if abs(x_new - x) < tol && abs(y_new - y) < tol
            x = x_new;
            y = y_new;
            break;
        end

        % Update x and y for the next iteration
        x = x_new;
        y = y_new;
        count = count + 1;
    end

    if count == max_iter
        warning('Maximum number of iterations reached without convergence.');
    end
end
