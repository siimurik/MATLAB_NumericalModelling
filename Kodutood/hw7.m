clc
clear

fprintf("Ülesanne 1.\n")
p = @(x) 4.*x.^5 - 2.*x.^4 + x.^3 + 7.*x - 9.0;

fprintf("a)\n")
x = [0, 6];
disp(p(x))

fprintf("b) Funktsiooni p(x) graafik\n")
figure(1)
fplot(p, [-1, 8])
title("p(x) = 4*x^5 - 2*x^4 + x^3 + 7*x - 9.")
grid on

fprintf("\nÜlesanne 2.\n")
% p(x) = x^2 − x − 12
coef = [1, -1, -12];
r = roots(coef);
disp(r)

fprintf("\nÜlesanne 3.\n")

roots = [-1, 3, 2, 5];
% Polynomial with specified roots or characteristic polynomial
p = poly(roots);
x = linspace(min(roots)-1, max(roots)+1, 1000);
% Polynomial evaluation
y = polyval(p, x);

figure(2);
plot(x, y);
hold on;

% Plot the roots
plot(roots, zeros(size(roots)), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

xlabel('x');
ylabel('p(x)');
title('Polünoom ja selle juured');
grid on;
hold off;

fprintf("\nÜlesanne 4.\n")

% Define the polynomials
p1 = [1 -2 1]; % Coefficients of p1(x) = x^2 - 2x + 1
p2 = [1 2];    % Coefficients of p2(x) = x + 2

% Convolution and polynomial multiplication
p = conv(p1, p2);
%{
MATLABi funktsioon conv() teostab konvolutsiooni, 
mis on matemaatiline operatsioon, mida kasutatakse 
sageli signaalitöötluses ja polünoomide korrutamisel. 
Kui kasutame conv() funktsiooni kahe vektori puhul, 
mis esindavad polünoomi koefitsiente, siis see on 
samaväärne nende polünoomide korrutamisega.
%}

% Display the resultant polynomial coefficients
disp('Resultant polynomial coefficients:');
disp(p);

% Define the range for x values for plotting
x = linspace(-3, 3, 1000);

% Evaluate the polynomial at each x value
y = polyval(p, x);

% Plot the polynomial
figure(3);
plot(x, y);
xlabel('x');
ylabel('p(x)');
title('Resultant Polynomial p(x)');
grid on;

fprintf("\nÜlesanne 5.\n")

x = [-1, 0.1, 0.5, 3, 4, 6.3,   7, 9, 14, 21];
y = [-2, 0.4, 0.7, 2, 4, 3.6, 3.8, 6, -1, 12];

% a) Fit a first-degree polynomial
p1 = polyfit(x, y, 1);

% b) Fit a second-degree polynomial
p2 = polyfit(x, y, 2);

% c) Fit a third-degree polynomial
p3 = polyfit(x, y, 3);

% d) Fit a sixth-degree polynomial
p6 = polyfit(x, y, 6);

% Define the range for x values for plotting
x_fit = linspace(min(x), max(x), 1000);

% Evaluate the polynomials at each x value
y1 = polyval(p1, x_fit);
y2 = polyval(p2, x_fit);
y3 = polyval(p3, x_fit);
y6 = polyval(p6, x_fit);

% Plot the data points and the fitted polynomials
figure(4);
plot(x, y, 'ko', 'MarkerFaceColor', 'k'); % Data points
hold on;
plot(x_fit, y1, 'r-', 'DisplayName', '1st Degree');
plot(x_fit, y2, 'g-', 'DisplayName', '2nd Degree');
plot(x_fit, y3, 'b-', 'DisplayName', '3rd Degree');
plot(x_fit, y6, 'm-', 'DisplayName', '6th Degree');
xlabel('x');
ylabel('y');
title('Polynomial Fits');
legend( "Location", "southwest");
grid on;
hold off;

% e) Find the polynomial value at x = 3.5 for 2nd and 6th degree polynomials
x_val = 3.5;
y_val_p2 = polyval(p2, x_val);
y_val_p6 = polyval(p6, x_val);

% Display the results
disp(['Value of 2nd degree polynomial at x = 3.5: ', num2str(y_val_p2)]);
disp(['Value of 6th degree polynomial at x = 3.5: ', num2str(y_val_p6)]);
