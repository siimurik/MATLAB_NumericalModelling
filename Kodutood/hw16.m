clc
clear

fprintf("Heuni meetod:\n")
f = @(x,y) exp(-y) + x;
h = 0.01;
x0 = 1.0;
xn = 4.5;
x = x0:h:xn;
y = zeros(1, length(x));
y(1) = 1.0; % algtingimus
for i = 2:length(x)
    y(i) = y(i-1) + h/2*( f(x(i-1), y(i-1)) + ...
            f(x(i), y(i-1) + h*f(x(i-1),y(i-1))));
end
fprintf("\ny'(%.1f) = %f\n", xn, y(end))
hold on
plot(x, y, "o-")

fprintf("Runge-Kutta meetod:\n")
steps_size = 0.05;
x_span = x0:steps_size:xn;
y0 = y(1);
[x,y] = ode45(f, x_span, y0);

fprintf("\ny'(%.1f) = %f\n", xn, y(end))
hold on
plot(x, y, "*-")