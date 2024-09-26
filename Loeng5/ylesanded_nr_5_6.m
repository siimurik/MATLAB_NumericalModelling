clc
clear

% Harjutus 5 ja 6, ülesanne 1
% vektrofunkstiooni defineerimine
fprintf("\nÜlesanne 1.")
F = @(x) [x(1)^3 - x(1)/x(2), sin(x(1)*cos(x(2)))]; % mittelineaarne võrrandisüsteem
F([2, 3])   % Kuna tegemist on vektrofunktsiooniga, siis sisend on samuti vektor
F([-3, -2])

fprintf("Ülesanne 2.\n")
F  = @(z) [z(1) - z(2)^2 + 1, -z(1)^2 + z(2)];   
x0 = [0, 0];
fsolve(F, x0)

% Kontroll
F(fsolve(F, x0))

% 2) lahendame hariliku.it. meetodiga
G = @(z) [-sqrt(z(2)), sqrt(z(1)+1)];
syms x y
dG = matlabFunction(jacobian([sqrt(y), sqrt(x+1)], [x, y]));
disp(dG)
norm(dG(-0.7, 0.5), 1) % < 1, sobib
norm(dG( 1.2, 1.4), 1) % < 1, sobib

z  = [-0.7, 0.5];
count  = 0;
max_iter = 1000;
% Iteratsiooniprotsess
while norm(F(z), 1) >= 1E-6 && count < max_iter
    z = G(z);
    count  = count + 1;
end
z
count

G = @(z) [sqrt(z(2)), sqrt(z(1)+1)];
z = [1.2, 1.4];
tol = 1E-6;
[~, ~] = HIM(F, G, z, tol, max_iter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minu meetod
figure(1)
fimplicit(@(x,y) x - y.^2 + 1)
hold on
fimplicit(@(x,y) -x.^2 + y)
legend("x - y^2 + 1", "-x^2 + y")
grid on

%g1 = @(x, y) (y.^2 -1);    % HAJUB
%g2 = @(x, y) (x.^2);       % HAJUB
g2 = @(x, y) sqrt(x+1);     % pos. sobib mõlema punkti jaoks
g1 = @(x, y) sqrt(y);       % neg. esimese ja pos. teise jaoks

g1_funcs = {@(x, y) -sqrt(y),   @(x, y) sqrt(y)};
g2_funcs = {@(x, y) sqrt(x+1),  @(x, y) sqrt(x+1)};

% Algväärtused
xy_vec = {[-0.7; 0.5], [1.2; 1.4]};
max_iter = 1000;
epsilon = 1E-6;

for i = 1:length(g1_funcs)
    g1 = g1_funcs{i};
    g2 = g2_funcs{i};
    xy = xy_vec{i};
    
    fprintf('Set %d: ', i);
    [~, ~, ~] = him(g1, g2, xy, epsilon, max_iter);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z, count] = HIM(F, G, x_vec, epsilon, max_iter)
    z = zeros(1, length(x_vec));
    for i = 1:length(x_vec)
        z(i) = x_vec(i);
    end
    count  = 0;
    % Harilik iteratsioonimeetod
    while norm(F(z), 1) >= epsilon && count < max_iter
        z = G(z);
        count  = count + 1;
    end
    if count == max_iter
        warning('Maximum number of iterations reached without convergence.');
    end
    %disp(z)
    %disp(count)
    fprintf('Iterations = %d, ', count)
    fprintf("Solution z = [%f, %f]\n", z(1), z(2))
end

function [x, y, count] = him(g1, g2, xy_vec, epsilon, max_iter)
    x = xy_vec(1);
    y = xy_vec(2);
    xvana = xy_vec(1) - 1;
    yvana = xy_vec(2) - 1;
    count = 0;
    while max(abs(x-xvana), abs(y-yvana)) >= epsilon && count < max_iter
        xvana = x;
        yvana = y;
        x = g1(x,y);
        y = g2(xvana,y);
        count = count + 1;
    end
    if count == max_iter
        warning('Maximum number of iterations reached without convergence.');
    end
    fprintf('Iterations = %d, Solution = [%.6f, %.6f]\n', count, x, y);
end