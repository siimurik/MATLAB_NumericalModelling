% Ülesanded 4, ülesanne 1
% dif. võrrandi numbriline lahendamine
fprintf("a) täpselt\n")
syms y(x);
difvor = diff(y, x) == y - 2.0*x;
lahend = matlabFunction(dsolve(difvor, y(0)==1))
figure(1)
fplot(lahend, [0, 2])
grid on
y_tapne = lahend(2)

f = @(x,y) y - 2.0.*x;
% y(0) = 1.0
[x, y] = ode45(f, [0, 2], 1.0)  % Runge-Kutta

hold on 
plot(x, y)

h = 0.4; % sammu pikkus
x = 0:h:2;

fprintf("b) Euleri meetod:\n")
y = zeros(1, length(x));
y(1) = 1.0; % algtingimus
for i = 2:length(x)
    y(i) = y(i-1) + h*f(x(i-1), y(i-1));
end
y(length(x))
hold on
plot(x, y, "o-")

fprintf("c) Heuni meetod:\n")
y = zeros(1, length(x));
y(1) = 1.0; % algtingimus
for i = 2:length(x)
    y(i) = y(i-1) + h/2*( f(x(i-1), y(i-1)) + ...
            f(x(i), y(i-1) + h*f(x(i-1),y(i-1))));
end
y(length(x))
hold on
plot(x, y, "*-")

fprintf("d) Keskpunkti meetod:\n");
y = zeros(1, length(x));
y(1) = 1.0; % algtingimus
for i = 2:length(x)
    y(i) = y(i-1) + h*f(x(i-1) + h/2, y(i-1) + h/2*f(x(i-1), y(i-1)));
end
y(length(x))
hold on
plot(x, y, ".--")

fprintf("e) Runge-Kutta meetod:\n");
y = zeros(1, length(x));
y(1) = 1.0; % algtingimus
c1    = 1.0; 
c2    = 1.0;
alpha = 0.45; 
beta  = 0.45;
for i = 2:length(x)
    y(i) = y(i-1) + c1*h*f(x(i-1), y(i-1)) + c2*h*f(x(i-1) + ...
                alpha*h, y(i-1) + beta*h*f(x(i-1), y(i-1)));
end
y(length(x))
hold on
plot(x, y, "-")