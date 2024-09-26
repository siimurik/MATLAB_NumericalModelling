clc
clear

fprintf("Ülesanne 1.\n")
f = @(x) 1.0./(1.0+x.^2);

a =  0.0;
b = 12.0;

fprintf("a)\n")
dx = (2-0)/(6-1)
x = 0:dx:2 % or x = linspace(0, 2, 6)
y = f(x)

fprintf("integraal vasakpoolse ristkülikvalemiga\n")
summa=0;
for i=1:5
    summa = summa + dx*y(i);
end
vasakp=summa

fprintf("integraal parempoolse ristkülikvalemiga\n")
summa=0;
for i=2:6
 summa = summa + dx*y(i);
end
paremp=summa


fprintf("e) täpne integraali\n")
integral(f, a, b)

%--------------------------------------------------------------------------
fprintf("Ülesanne 2.\n")
x = [0.00 0.25 0.50 0.75 1.00];
f = [4.32 4.36 4.58 5.79 6.14];
plot(x, f, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r')

fprintf("a) integraal trapets ja Simpsoni valemiga\n")
trapets = trapz(x, f)

% a) Simpsoni meetodi jaoks lähendame funktsiooni kuupsplainiga S^(3,2)(x)
x_vec = [0:0.001:1]; % suurendame sõlmede arvu, et splainifunktsioon

% oleks joonisel sujuvam
y_cubeVec = interp1(x, f, x_vec, 'cubic'); % vector of f(xx) values 
hold on
plot(x_vec, y_cubeVec)
grid on

% esitame kuupsplaini funktsioonina
func3 = @(z) interp1(x, f, z, 'cubic');
%func2 = @(z) interp1(x, f, z,'cubic')


% c) leiame integraali Simpsoni valemiga
fprintf("Kuuppolünoomi tulemus:\n")
integral(func3, 0, 1)    % built-in

% True Simpson. Janno lk. 186
simpson = 0;
a = 0.0;
b = 1.0;
n = 1000;
h = (b-a)/n;
for i = 1:n/2
    simpson = simpson + y_cubeVec(2*i-1) + 4*y_cubeVec(2*i) + y_cubeVec(2*i+1);
end
simpson = h/3*simpson


fprintf("Ruutpolünoomi tulemus:\n")
coefs = polyfit(x, f, 2);
func2 = @(x) coefs(1).*x.^2 + coefs(2).*x + coefs(3);
y_squarVec = func2(x_vec);

integral = simpsonSolve(y_squarVec, 0, 1)

fprintf("Manual kuuppolünoomi tulemus:\n")
coefs = polyfit(x, f, 3);
f3 = @(x) coefs(1).*x.^3 + coefs(2).*x.^2 + coefs(3).*x + coefs(4);
cubevec = f3(x_vec);

integral = simpsonSolve(cubevec, 0, 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function integral = simpsonSolve(funcVec, a, b)
    n = length(funcVec)-1;
    integral = 0;
    if mod(n,2) == 0
        h = (b-a)/n;
    else
        n = n-1;
        h = (b-a)/n;
    end
    for i = 1:n/2
        integral = integral + funcVec(2*i-1) + 4*funcVec(2*i) + funcVec(2*i+1);
    end
    integral = h/3*integral;
end
