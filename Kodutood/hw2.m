clc
clear
% Iseseisev töö nr. 2

disp("Ülesanne 1.")
f = @(x) 3.*x.^4 - x.^3 + 5.*x.^2 + 7;
x1 = -1.0;
x2 =  4.0;
fprintf("f(%d) = %f\n", x1, f(x1));
fprintf("f(%d) = %f\n", x2, f(x2));

x = linspace(-2, 5);
y = f(x);

figure(1);
plot(x, y, 'b-', 'LineWidth', 1);
hold on;

% Punktid (-1, f(-1)) ja (4, f(4))
plot(-1, f(-1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(4, f(4), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
title('Funktsiooni y = 3x^4 - x^3 + 5x^2 + 7 graafik')
xlabel('x')
ylabel('y')
grid on
hold off

disp("Ülesanne 2.")
f = @(x) x - cos(x);
df = @(x) 1 + sin(x);

x = linspace(0, 2*pi);

y = f(x);
dy = df(x);

figure(2);
plot(x, y, 'b-') % Funktsiooni graafik
hold on
plot(x, dy, 'r--') % Tuletise graafik

title('Funktsiooni y = x - cos(x) ja selle tuletise graafikud');
xlabel('x')
ylabel('y')
legend('y = x - cos(x)', 'y'' = 1 + sin(x)')
grid on
hold off

disp("Ülesanne 3.")
clear
figure(3)
%ezplot(@(x,y) x.^2.*y - 4.*y.^2 + 3)    % töötab sest funktsioon on ilmutamata
fimplicit(@(x,y) x.^2.*y - 4.*y.^2 + 3)
grid on
title('Funktsiooni x^2 y - 4 y^2 + 3 = 0 graafik')
xlabel('x')
ylabel('y')
hold off

figure(4)
x = [-1, 2, 3, 5];
y = [4, -6, 0, 3];
plot(x, y, 'bo-', 'MarkerFaceColor', [0 0.7 0.7])
grid on
hold off

figure(5)
x = @(t) sin(t);
y = @(t) cos(t);
title('Funktsioonide x = sin(t) ja y = cos(t) graafik')
fplot(x, y)
grid on
hold off

figure(6)
y = @(x) x.^4 - 3.*x.^2 + x + 1;
fplot(y)
title('Funktsiooni y = x^4 - 3*x^2 + x + 1 graafik')
grid on
hold off

figure(7)
%ezplot(@(x,y) y.^2 - 3.*x - 2)
fimplicit(@(x,y) y.^2 - 3.*x - 2)
grid on
title('Funktsiooni y^2 - 3*x - 2 = 0 graafik')
hold off

disp("Ülesanne 4.")
%figure(8)
y = @(x) 2.*x.^3 + 4.*x.^2 + 6.*x -7;
fplot(y, [-3, 4])
grid on
hold off

disp("Ülesanne 5.")
x = linspace(-10, 10, 400);
y = linspace(-10, 10, 400);
[X, Y] = meshgrid(x, y);

F1 = 3*X.^2 + Y.^2 - 8;
F2 = (X-4).^2 + (Y-3).^2 - 16;

figure(9);
contour(X, Y, F1, 'r')
hold on
contour(X, Y, F2, 'b')
hold off

legend('3x^2 + y^2 = 8', '(x-4)^2 + (y-3)^2 = 16')
title('Funktsioonide graafikud')
xlabel('x')
ylabel('y')
grid on
hold off

disp("Ülesanne 6.")
log3 = @(x) log(x) / log(3); % https://en.wikipedia.org/wiki/Logarithm#Change_of_base

disp("a)")
z = @(x,y) log3(x) + 4*exp(cos(5.*x+sqrt(y)))

disp("b)")
x = 7;
y = 1;
fprintf("x = %f\n", x);
fprintf("y = %f\n", y);
fprintf("z(x,y) = %f\n\n", z(x,y));

disp("c)")
fprintf("z(%d,%d) = %f\n\n", 6, 2, z(6,2));

disp("d)")
figure(10)
x = linspace(0, 5);
y = linspace(0, 5);
[X, Y] = meshgrid(x,y);
Z = z(X, Y);
surf(X, Y, Z)
grid on
hold off

disp("Ülesanne 7.")
alpha = 4.;
beta  = 5.;
gamma = 1.;

a = sqrt((exp(alpha) + 2*beta) / ((alpha - beta)^4 + beta^alpha));
b = (alpha - beta) / (alpha + beta) * gamma;

f = @(x) (x.^2 - 2.*a.*b.*x) ./ (x.^2 + 2.*a.*b.*x) .* cos(x);
figure(11)
fplot(f, [-10, 10]);
xlabel('x');
ylabel('f(x)');
title('Plot of the function f(x)');
grid on;
hold off

disp("Ülesanne 8.")
% Funktsioon x^3 + 5*x + 4 = 0
% aka 1.0*x^3 + 0*x^2 + 5*x + 4 = 0
coef = [1 0 5 4];
r = roots(coef)
figure(12)                          % NOTE:
%ezplot(@(x) x.^3 + 5.*x + 4)       % töötab. // Lisafunktsioon: lisab ka 'title' vastava funktsiooni jaoks.
fplot(@(x) x.^3 + 5.*x + 4)         % töötab
%fimplicit(@(x) x.^3 + 5.*x + 4)    % ei tööta... miks???
grid on
fprintf("Graafikult võimalik näha nullkohta [%f, %f].\n", real(r(3)), 0.0);
hold on 
plot(real(r(3)), 0.0, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
fprintf("\n");
hold off

disp("Ülesanne 9.")
% Funktsioon 3*x^6 + 2*x^5 + 12*x^4 - 7*x^2 - x - 16 = 0
% aka 3*x^6 + 2*x^5 + 12*x^4 + 0*x^3 - 7*x^2 - x - 16 = 0
coef = [3 2 12 0 -7 -1 -16];
r = roots(coef)
f = @(x) 3.*x.^6+2.*x.^5+12.*x.^4-7.*x.^2-x-16;

figure(13)
fplot(f, [-2, 2], 'b-');
ylim([-60, 100]);
grid on

fprintf("Nullkohad antud funktsiooni jaoks:\n");
fprintf("3*x^6 + 2*x^5 + 12*x^4 - 7*x^2 - x - 16 = 0\n");

% Reaalarvulised nullkohad
x_zeros = [];
for i = 1:length(r)
    if imag(r(i)) == 0.0
        fprintf("[%f, %f]\n", real(r(i)), 0.0)
        x_zeros = [x_zeros, real(r(i))];
    end
end
% Iga reaalarvulise nullkoha lisamine graafikule
hold on
for i = 1:length(x_zeros)
    plot(x_zeros(i), 0.0, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
end
hold off
fprintf("\n");


disp("Ülesanne 10.")
% Define matrix A and vector b
A = [67 12 16; -54 71 43; -82 35 88];
b = [31; 15; -14];
%linsolve(A, b)  % Too easy

% The supercool custom BiCGSTAB solver:
% Source: https://ctk.math.ncsu.edu/matlab_roots.html
% https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
% Step 1:
% Define initial guess x0 (zero vector)
x0 = zeros(3, 1);

% Step 2:
% A matrix-vector product routine
%   atv must return Ax when x is input
%	the format for atv is
%   function ax = atv(x)
atv = @(x) A*x;

% Step 3:
% Set the solver parameters
params = [1e-6, 100];  % params(1): tolerance, params(2): max iterations

% Step 4:
% Call the custom Bi-CGSTAB solver
[x, error, total_iters] = BiCGSTAB(x0, b, atv, params);

% Built-in version. Slower.
%[x, flag, error, total_iters] = bicgstab(A, b, params(1), params(2)); 

% Display the results
disp('Solution x:');
disp(x);

disp('Error at each iteration:');
fprintf("%e\n", error);
%disp(error);

disp('Total number of iterations:');
disp(total_iters);

disp("Ülesanne 11.")
A = [4.  5. -3. 1.; 
     2.  0. -1. 7.; 
     3.  4. -6. 2.;
     7. -1. -7. 3.];
b = [2.; -5.; 8.; 9.];
%tic; x = A \ b; toc;
x0       = zeros(4, 1);
atv      = @(x) A*x;
tol      = 1.0E-6;
max_iter = 100; 
params   = [tol, max_iter];
%[x, error, total_iters] = BiCGSTAB(x0, b, atv, params);

% Transpose-Free Quasi-Minimal Residual method
% Source: https://ctk.math.ncsu.edu/matlab_roots.html
% Explanation: https://mathworld.wolfram.com/Quasi-MinimalResidualMethod.html
[x, error, total_iters] = TFQMR(x0, b, atv, params);

fprintf('Solution x:\n');
disp(x)
fprintf('Error at each iteration:\n');
fprintf("%e\n", error);
fprintf('Total number of iterations: %d\n', total_iters);

