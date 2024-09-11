% Funktsioon graafikud, näide
% funktsioon tabeli kujul
x = [0 1 2 3];
y = [0 1 4 9];
plot(x, y,'go-')

%b) funktsioon valemi kujul, y = f(x)
hold on
fplot(@(x) x.^2)

%c) funktsioon ilmutatmata kujul, st
% f(x,y) = 0
hold on
fimplicit(@(x,y) x.^2 + y.^2 - 1) 
axis equal % vastasel juhul näeb välja nagu ellips