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
hold off

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


fprintf("Ülesanne 3.\n")
figure(2)
fimplicit(@(x,y) x.^2 + 2.*y - 8)
hold on
fimplicit(@(x,y) x + y.^2 - y - 3)
legend("x.^2 + 2.*y - 8", "x + y.^2 - y - 3")
grid on
hold off
% alglähendid graafiku alusel:
%   * (-3.8, -2.1)
%   * (-1.6,  2.7)
%   * ( 2.2,  1.4)
%   * ( 2.8, -0.1)

F = @(z) [ z(1).^2 + 2.*z(2) - 8,  z(1) + z(2).^2 - z(2) - 3];
x0 = [0, 0];
fsolve(F, x0)

% lahendame hariliku.it. meetodiga
G = @(z) [ -sqrt( 8.0 - 2.*z(2) ), -sqrt(z(2) - z(1) + 3.0)];
syms x y
dG = matlabFunction(jacobian([ sqrt( 8.0 - 2.*y ),  sqrt(y - x + 3.0)], [x, y]));
disp(dG)
alglahend = {[-3.8, -2.1], [-1.6,  2.7], [ 2.2,  1.4], [ 2.8, -0.1]};
for i = 1:length(alglahend)
    z = alglahend{i};
    norm(dG(z(1), z(2)), 1) % < 1, sobib
end
% viimane ei sobi
% Leiame need mis sobisid
z        = alglahend{1};
count    = 0;
max_iter = 1000;
% Iteratsiooniprotsess
while norm(F(z), 1) >= 1E-6 && count < max_iter
    z = G(z);
    count  = count + 1;
end
z
count

% Natuke mugavam viis 
G_funcs = {@(z)[ -sqrt( 8.0 - 2.*z(2) ), -sqrt(z(2) - z(1) + 3.0)], ...
           @(z)[ -sqrt( 8.0 - 2.*z(2) ),  sqrt(z(2) - z(1) + 3.0)], ...
           @(z)[  sqrt( 8.0 - 2.*z(2) ),  sqrt(z(2) - z(1) + 3.0)]};
%           @(z)[  sqrt( 8.0 - 2.*z(2) ), -sqrt(z(2) - z(1) + 3.0)]}; % Pole vaja lisada sest nagunii ei koondu.

tol = 1.0E-6;
[~, ~] = HIM_Vec(F, G_funcs, alglahend, tol, max_iter);

% 4. alglähendiga teeme nata sohki kuna 1-norm dG tuli 1
G = @(z) [ sqrt( 8.0 - 2.*z(2) ), z(1) + z(2).^2 - 3]; % defineerime uue G funktsiooni
z        = alglahend{4};
count    = 0;
max_iter = 1000;
% Iteratsiooniprotsess
while norm(F(z), 1) >= 1E-6 && count < max_iter
    z = G(z);
    count  = count + 1;
end
z
count

%-------------------------------------------------------------
% Newtoni meetodiga
F  = @(z) [ z(1).^2 + 2.*z(2) - 8;
            z(1) + z(2).^2 - z(2) - 3];
dF = @(z) [ 2.*z(1),    2.0; 
            1.0 ,       2.*z(2)-1.0];

z = [ 2.1 ; 
      1.6 ];
count = 0;
maxit = 1000;
while norm(F(z), 1) >= 1E-6
    z = z - dF(z)^(-1) * F(z); % (2x1) - (2x2)*(2x1)
end
z
count

initial_values = {[-3.8; -2.1], [-1.6;  2.7], [2.2;  1.4], [ 2.8; -0.1]};
tol = 1E-6;
disp("Newtoni meetod:")
[~, ~] = newton_method(F, dF, initial_values, tol, max_iter);
disp("Modifitseeritud Newtoni meetod:")
[~, ~] = modified_newton(F, dF, initial_values, tol, max_iter);

fprintf("Ülesanne 4.\n")
figure(3)
fimplicit(@(x,y) 20.*x.^2 - 1.0 + 2.*x.^3 - 4.*y.^3)
hold on
fimplicit(@(x,y) 2.*x.^3 - 3.*y.^2 - 10.*y + 5.0)
legend("20.*x.^2 - 1.0 + 2.*x.^3 + 4.*y.^3", "2.*x^3 - 3.*y.^2 - 10.*y + 5.0")
grid on
hold off
% alglähendid graafiku alusel:
    % * (-0.2,  0.4)
    % * ( 0.2,  0.4)
    % * ( 3.6,  4.5)

syms x y
dG = matlabFunction(jacobian([ ...
%                                sqrt(1/20*(1.0 - 2.*x.^3 + 4.*y.^3)),  ...  % 1. ja 2. punkt
%                                1/10*(2*x^3 - 3*y^2 + 5)]           ,       % 1. ja 2. punkt
%                                nthroot( 1/2*(1.0 - 20*x^2 + 4*y^3), 3), ...
%                                    1/10*(2*x^3 - 3*y^2 + 5)
%                                sqrt(1/3*(2*x^3 - 10*y + 5)) ...       % 
                                nthroot( 1./4*(20.*x.^2 - 1.0 + 2.*x.^3), 3) ...    % 3. punkt
                                nthroot( 1/2*(3.*y.^2 + 10.*y - 5.0), 3), ...       % 3. punkt
                                                                  ],   [x, y]));


disp(dG)
alglahend = {[-0.2,  0.4], [0.2,  0.4], [3.6,  4.5]};
for i = 1:length(alglahend)
    z = alglahend{i};
    norm(dG(z(1), z(2)), 1) % < 1, sobib
end

G_funcs = { @(z)[ -sqrt(1./20*(1.0 - 2.*z(1).^3 + 4.*z(2).^3)), 1./10*(2*z(1).^3 - 3.*z(2).^2 + 5) ], ...
            @(z)[ +sqrt(1./20*(1.0 - 2.*z(1).^3 + 4.*z(2).^3)), 1./10*(2*z(1).^3 - 3.*z(2).^2 + 5) ], ...
            @(z)[ nthroot(1/2*(3.*z(2).^2 + 10.*z(2) - 5.0), 3), nthroot( 1./4*(20.*z(1).^2 - 1.0 + 2.*z(1).^3), 3)]};

alglahend = {[-0.2,  0.4], [0.2,  0.4], [3.6,  4.5]};
tol = 1.0E-6;
max_iter = 100000;
[~, ~] = HIM_Vec(F, G_funcs, alglahend, tol, max_iter);
% Praegu funktsioonid ei tundu koonduvat. Üks trikk on funkstioonide mõlemaid pooli liita x-ga läbi ehk
% x + 20.*x.^2 = 1.0 - 2.*x.^3 + 4.*y.^3 + x

%--------------------------------------------------------------
% Newtoni meetodiga
clear
fprintf("\nLahendite leidmine Newtoni ja mod. Newtoni meetodiga:")
F  = @(z) [     20.*z(1).^2 - 1.0 + 2.*z(1).^3 - 4.*z(2).^3;
                2.*z(1).^3 - 3.*z(2).^2 - 10.*z(2) + 5.0];
syms x y
fprintf("Vektorfunktsiooni F tuletis:\n")
jacobian([20.*x.^2 - 1.0 + 2.*x.^3 - 4.*y.^3 2.*x.^3 - 3.*y.^2 - 10.*y + 5.0], [x,y])
dF = @(z)[  6.*z(1).^2 + 40.*z(1), -12.*z(2).^2;
            6.*z(1).^2,           -6.*z(2) - 10];

alglahend = {[-0.2,  0.4], [0.2,  0.4], [3.6,  4.5]};
tol = 1E-6;
max_iter = 1000;
disp("Newtoni meetod:")
[~, ~] = newton_method(F, dF, alglahend, tol, max_iter);
disp("Modifitseeritud Newtoni meetod:")
[~, ~] = modified_newton(F, dF, alglahend, tol, max_iter);

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

function [z, count] = HIM_Vec(F, G_funcs, xy_vec, epsilon, max_iter)
    for i = 1:length(G_funcs)
        z = xy_vec{i};
        G = G_funcs{i};
        count  = 0;
        % Harilik iteratsioonimeetod
        while norm(F(z), 1) >= epsilon && count < max_iter
            z = G(z);
            count  = count + 1;
        end
        if count == max_iter
            warning('Maximum number of iterations reached without convergence.');
        end
        fprintf('Iterations = %d, ', count)
        fprintf("Solution z = [%f, %f]\n", z(1), z(2))
    end
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

function [z, count] = newton_method(F, dF, initial_values, epsilon, max_iter)
    % F - funktsioonide veeruvektor, mille nullkohta otsime
    % dF - jakobiaani maatriks
    % initial_values - alglähendide paarid / ristumispunktid
    % epsilon - täpsus
    % max_iter - maksimaalne iteratsioonide arv

    % Käime läbi iga arvupaari
    for i = 1:length(initial_values)
        z = initial_values{i};
        count = 0;
          
        while norm(F(z), 1) > epsilon && count < max_iter
            z0 = z;
            z  = z0 - dF(z0)^(-1) * F(z0); % Newtoni meetodi eeskiri
            count = count + 1;
        end
          
        fprintf('Set %d: Iterations = %d, Solution = [%.6f, %.6f]\n', i, count, z(1), z(2));
    end
end

function [z, count] = modified_newton(F, dF, initial_values, epsilon, max_iter)
    % F - funktsioonide veeruvektor, mille nullkohta otsime
    % dF - jakobiaani maatriks
    % initial_values - alglähendide paarid / ristumispunktid
    % epsilon - täpsus
    % max_iter - maksimaalne iteratsioonide arv

    % Käime läbi iga arvupaari
    for i = 1:length(initial_values)
        z = initial_values{i};
        z0 = z;
        count = 0;
          
        while norm(F(z), 1) > epsilon && count < max_iter
            z = z - dF(z0)^(-1) * F(z); % modifitseeritud Newtoni meetodi eeskiri
            count = count + 1;
        end
          
        fprintf('Set %d: Iterations = %d, Solution = [%.6f, %.6f]\n', i, count, z(1), z(2));
    end
end