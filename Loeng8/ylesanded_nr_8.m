clc
clear

fprintf("\n Ülesanne 6.\n")
syms x y
f = @(x,y) cos(x*y) + 3*exp(sin(7*x)) - nthroot( 2*x^3*y^4+4, 5);
diff(f, x)
diff(f, y)

fprintf("\nÜlesanne 3.\n")
A = [-5 3; 2 -1];
B = [4 -2; 7 0];

fprintf("\na) C maatriks, mis on A ja B elementite korrutis")
C = A.*B

fprintf("\nb) B maatrikis teise rea esimene element")
B(2, 1)

fprintf("\nÜlesanne 3.\n")
clear
y = @(x) x^3 - 2*x^2 + 7;
y(5)
% Või
x = 5;
y(x)

x = [0, 1, 4, 3];
y = @(x) x.^3 - 2.*x.^2 + 7;
y(x)

fprintf("\nÜlesanne 1.\n")
fprintf("\na) graafik MATLABis")
figure(1)
f = @(x) x.^5 - 2.*x.^3 + 2.*x - 2;
fplot(f)
grid on

fprintf("\nb) graafik MATLABis")
syms x
g = @(x) nthroot(2*x^3 - 2*x + 2, 5);
%g = @(x) nthroot(0.5*(x^5 + 2*x -2), 3);
dg = matlabFunction(diff(g, x))
abs(dg(1.1)) 

fprintf("\nc) x1 ja x2")
clear
x=1.1;
g = @(x) nthroot(2*x^3 - 2*x + 2, 5);
for i = 1:2
    fprintf("x(%d) = g(%f) = %f", i, x, g(x));
    x = g(x)
end

fprintf("\nd) Võrrandi lahend")
f = @(x) x.^5 - 2.*x.^3 + 2.*x - 2;
fsolve(f, 1.1)

fprintf("\nÜlesanne 5.\n")
figure(2)
fimplicit(@(x,y) x.^2 - 2.*y + 4)
hold on
fimplicit(@(x,y) y.^2 - x.^4 + 4.*x - 1)
title("Võrrandisüsteemi funktsioonid")
xlabel("x-telg")
ylabel("y-telg")
legend("x^2 - 2*y + 4 = 0", "y^2 - x^4 + 4*x - 1 = 0")