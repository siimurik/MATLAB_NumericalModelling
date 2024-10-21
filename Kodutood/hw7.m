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

fprintf("\nÜlesanne 6.\n")
x = [-1.0  0.0 1.0  2.0];
y = [ 3.0 -4.0 5.0 -6.0];
lagrp = lagrange_poly(x, y);
fprintf("Lagrange's interpolation polynomial Φ(x) =")
disp(lagrp)

fprintf("f(0.5) = %.4f\n", lagrp(0.5))

fprintf("\nÜlesanne 7.\n")
x_aasta = [2014 2015 2016 2017 2018 2019 2020 2021];
y_oC    = [22.2 19.5 18.0 17.9 25.2 18.6 18.7 29.3];
figure(5)
plot(x_aasta, y_oC, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
grid on
hold on

% 1st-degree polynomial fit
fprintf("Lähendamine 1. astme polünoomiga:\n")
p1 = polyfit(x_aasta, y_oC, 1);   % coefficients of the polynomial
disp(p1)
xdata = linspace(min(x_aasta)-1, max(x_aasta)+1, 1000);
y1 = polyval(p1, xdata);
plot(xdata, y1)

% Centering and scaling x_aasta to avoid ill-conditioning for higher-degree polynomial
fprintf("Lähendamine 4. astme polünoomiga (centered and scaled):\n")
x_mean = mean(x_aasta);
x_scale = std(x_aasta);

% Scale and center x-values
x_scaled = (x_aasta - x_mean) / x_scale;

% Perform the polynomial fit on the scaled x-values
p4 = polyfit(x_scaled, y_oC, 4);
disp(p4)

% Generate the scaled x-values for plotting
xdata_scaled = (xdata - x_mean) / x_scale;
y4 = polyval(p4, xdata_scaled);
plot(xdata, y4)

% Add labels
xlabel('Year');
ylabel('Temperature (°C)');
legend('Data points', '1st-degree fit', '4th-degree fit');


%-----------------------------------------------------------------
function lagrp = lagrange_poly(x, y)
    % Validate input lengths
    if length(x) ~= length(y)
        error('Input vectors x and y must have the same length.');
    end
    
    % Get the number of points
    n = length(x);
    
    % Define symbolic variable
    syms X;
    
    % Initialize the Lagrange polynomial
    lagrp_sym = 0;
    
    % Construct the Lagrange polynomial
    for i = 1:n
        % Compute the Lagrange basis polynomial L_i(X)
        Li = 1;
        for j = 1:n
            if i ~= j
                Li = Li * (X - x(j)) / (x(i) - x(j));
            end
        end
        
        % Add the contribution of the current basis polynomial to the total
        lagrp_sym = lagrp_sym + Li * y(i);
    end
    
    % Simplify the resulting polynomial
    lagrp_sym = expand(lagrp_sym);
    
    % Convert to a MATLAB function for numerical evaluation
    lagrp = matlabFunction(lagrp_sym);
end
