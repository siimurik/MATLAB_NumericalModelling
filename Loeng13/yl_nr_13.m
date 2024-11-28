clc
clear

syms y(x)
difvor = 2.0 * x - sin(y)*diff(y,x) == 0
yldlahend = matlabFunction(dsolve(difvor))

C = [1 3 6];
figure(1)
fplot(yldlahend(C,x))
grid on

difvor3 = 2.0*x*y + (y-1)/exp(x)*diff(y,x) == 0
yldlahend3 = (dsolve(difvor3))

difvor4 = diff(y,x) - cot(x)*y - 2.0*x*sin(x) == 0
yldlahend4 = (dsolve(difvor4))

difvor5 = diff(y,x,2) - 2.0*diff(y,x) == 0
yldlahend5 = (dsolve(difvor5))

difvor6 = diff(y,x,2) + 4.0*diff(y,x) + 4.0*y == 0
yldlahend6 = (dsolve(difvor6))
erilahend6 = dsolve(difvor6, [y(0) == 1, subs(diff(y,x), x, 0) == 0])