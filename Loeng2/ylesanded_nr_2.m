clc
clear

% Ülesanne 1
figure(1)
fplot(@(x) x.^5 - x.^2 + 2, [-2 2])
title('Funkstiooni graafik')
legend('y=x^5 - x^2 + 2') % võib ka tühjaks jätta ja töötab
xlabel('x-telg')
ylabel('y-telg')
grid on
hold off

% Ülesanne 2
y1 = @(x)2.*x + 4;
y2 = @(x)x.^2 - 2.* + 1;
y3 = @(x) 1./x;
figure(2)
hold on
grid on
%fplot(y1, [0.01 10])
%fplot(y2, [0.01 10])
%fplot(y3, [0.01 10])
fplot(@(x) [2.*x + 4, x.^2 - 2.* + 1, 1./x], [0.01 10])
title('Funkstiooni graafik')
legend('y=2*x + 4', 'y=x^2 - 2* + 1', 'y=1/x')
xlabel('x-telg')
ylabel('y-telg')
hold off

% Ülesanne 3
t = 10.0;
y = @(x) 5.*sqrt(x.*t) + t*sin(x);
figure(3)
fplot(y, [0 10])
hold off

% Ülesanne 4
figure(4)
fimplicit(@(x,y) (x-1).^2 + (y+2).^2 - 25, [-6 8 -8 4])
xlabel('x-telg')
ylabel('y-telg')
grid on
legend('(x-1)^2 + (y+2)^2 = 25')
axis equal
hold off

% Ülesanne 5
figure(5)
hold on
fimplicit(@(x,y) x.^2 + y.^2 - 4, [-3 3 -4 4])
fimplicit(@(x,y) 10.*x.^2 + y.^2 - 9, [-3 3 -4 4])
legend('x^2 + y^2 = 4', '10*x^2 + y^2 - 9 = 0')
grid on
axis equal
hold off

% Ülesanne 6
A = [1  4  7
     8  2  0
     5 -1  1];

B = [6; 10; 3];
d = [1:1:10]';

F = transpose(A);

% Ax = F -> x = F*A^-1
x = A^(-1)*F
%x = linsolve(A, F)
%A \ F

A(1,2) + F(2,1) + d(10)

% Ülesanne 7
A = [4.  2.  0.
     2.  3.  1.
     0.  1.  2.5];

b = [2; 5; 6];

tic
A^(-1)*b
toc

tic
A \ b
toc

tic
linsolve(A, b)
toc