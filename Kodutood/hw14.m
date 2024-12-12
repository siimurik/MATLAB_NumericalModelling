clc
clear

fprintf("Ülesanne 1.\n")
%{
 dy         
---- = y - 2x
 dx         

 y(0) = 1
%}
syms y(x)
fprintf("a) Täpne lahend:\n")

func = diff(y,x,1) == y - 2.0*x;
genSol = dsolve(func);
fprintf("Üldlahend:\n"); disp(genSol)

specSol = dsolve(func, y(0) == 1); 
fprintf("Erilahend:\n"); disp(specSol)

y_eri = matlabFunction(specSol);
fprintf("f(2.0) = %f \n\n", y_eri(2.0))

fprintf("b)/c) Euleri meetodiga käsitsi ja kasutades MATLABi.\nSammu pikkus olgu h = 0.4 / h = 0.1;\n")

h = 0.1;    % <--- change out the value here
% y(0) = 1.0
x0 = 0.0;
y0 = 1.0;
xn = 2.0;

[x, y, yeuler] = euler_method(h, x0, xn, y0);

disp('x väärtused:');
disp(x);
disp('y väärtused:');
disp(y);
disp('Euler meetodi tulemus, kui dy/dx(2.0):');
disp(yeuler);

fprintf("Ülesanne 2.\n")
%{
 dy         
---- = y + 1
 dx

 y(0) = 1
%}
syms y(x)
fprintf("a) Täpne lahend:\n")

func = diff(y,x,1) == y + 1.0;
genSol = dsolve(func);
fprintf("Üldlahend:\n"); disp(genSol)

specSol = dsolve(func, y(0) == 1); 
fprintf("Erilahend:\n"); disp(specSol)

y_eri = matlabFunction(specSol);
fprintf("f(3.0) = %f \n\n", y_eri(3.0))

fprintf("b) Leida diferentsiaalvõrrandi lahendi väärtus, kui x = 3 Euleri meetodiga MATLABis, kasutades erinevaid sammu pikkuseid\n")

h = 0.0001;    % <--- change out the value here
% y(0) = 1.0
x0 = 0.0;
y0 = 1.0;
xn = 3.0;

[x, y, yeuler] = euler_method(h, x0, xn, y0);

disp("Sammu pikkus:")
disp(h)
% disp('x väärtused:');
% disp(x);
% disp('y väärtused:');
% disp(y);
disp('Euler meetodi tulemus, kui dy/dx(3.0):');
disp(yeuler);


%----------------------------------------------------------------------------------------------------
function [x, y, yeuler] = euler_method(h, x0, xn, y0)
    % Euleri meetodiga lahendamine
    n = (xn - x0) / h; % osalõikude arv
    x = zeros(1, n+1); % x väärtuste vektor
    y = zeros(1, n+1); % y väärtuste vektor
    x(1) = x0; % algtingimuse x-i väärtus
    y(1) = y0; % algtingimuse y-i väärtus
    f = @(x, y) y - 2.0 * x; % funktsioon f(x, y)
    
    for i = 2:(n+1)
        x(i) = x(i-1) + h;
        y(i) = y(i-1) + h * f(x(i-1), y(i-1));
    end
    
    yeuler = y(n+1); % lõppväärtus y
end