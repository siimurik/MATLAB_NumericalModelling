clc
clear

fprintf("Ülesanne 1.\n")
%{
   dy
x --- = y(y-1)
   dx
%}
syms y(x)
vorrand = x*diff(y, x, 1) == y*(y-1);
yldlahend = matlabFunction(dsolve(vorrand));
fprintf("Üldlahend:\n")
disp(yldlahend)

% b) esitame moned diferentsiaalvorrandi lahendid joonisel
y1 = yldlahend(1, x);
y2 = yldlahend(2, x);
y3 = yldlahend(3, x);
y4 = yldlahend(4, x);

figure(1)
fplot(y1(1), 'r')
hold on
fplot(y2(1), 'g')
fplot(y3(1), 'b')
fplot(y4(1), 'm')
legend("Location","southwest")
grid on
hold off
% c) erilahend algtingimusel y(1)=1
erilahend=dsolve(vorrand, y(1)==1);
fprintf("Erilahend:\n"); disp(erilahend)


fprintf("Ülesanne 2.\n")
%{
          dy
y*(x^2-1)--- = xy^2 + x
          dx
%}
syms y(x) % Define y as a symbolic function of x
func = y*(x^2-1)*diff(y,x,1) == x*y^2 + x; % Differential equation

% Solve the general solution symbolically
genSol = dsolve(func);
fprintf("Üldlahend:\n"); disp(genSol)

genSolFunc = matlabFunction(genSol);
y1 = genSolFunc(1, x);
y2 = genSolFunc(2, x);
y3 = genSolFunc(3, x);
y4 = genSolFunc(4, x);
figure(2)
fplot(y1(1), 'r')
hold on
fplot(y2(1), 'g')
fplot(y3(1), 'b')
fplot(y4(1), 'm')
xlim([-5, 5])
ylim([-200, 50])
legend("y1","y2","y3","y4","Location","south")
grid on
hold off

% Specify the initial condition symbolically
specSol = dsolve(func, y(2) == -1); 
fprintf("Erilahend:\n"); disp(specSol)

% Plot the specific solution (only real parts if necessary)
figure(3)
fplot(real(specSol), [-10, 10]) % Plot real part of the solution
grid on
ylim([-5, 3])

fprintf("Ülesanne 3.\n")
%{
xy dy + (x^2-1)dx = 0
or
   dy       (x^2 - 1)
y ---- = - -----------
   dx           x
%}
func = y*diff(y,x,1) == - (x^2-1)/x;
genSol = dsolve(func);
fprintf("Üldlahend:\n"); disp(genSol)

specSol = dsolve(func, y(1) == 0); 
fprintf("Erilahend:\n"); disp(specSol)

fprintf("Ülesanne 4.\n")
%{
 dy         2
---- - y = x
 dx         
%}
func = diff(y,x,1) - y == x^2;
genSol = dsolve(func);
fprintf("Üldlahend:\n"); disp(genSol)

fprintf("Ülesanne 5.\n")
%{
 dy         
---- = 2(2x - y)
 dx         
%}
func = diff(y,x,1) == 2*(2*x - y);
genSol = dsolve(func);
fprintf("Üldlahend:\n"); disp(genSol)

specSol = dsolve(func, y(0) == -1); 
fprintf("Erilahend:\n"); disp(specSol)
