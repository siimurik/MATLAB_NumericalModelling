clc
clear

fprintf("Ülesanne 1.\n");
fprintf("a)\n");
fa1 = @(x)  x.^2;   % x <= 0
fa2 = @(x) -x.^2;   % 0 < x <= 1
fa3 = @(x) 1-2.*x;  % x > 1

fprintf("\nKontrollime, kas esitab splaini (kui saame 0, siis ei):\n")
disp(fa1(0.0) == fa2(0.0))
disp(fa2(1.0) == fa3(1.0))
fprintf("Kõik on tõesed seega edaspidi leiame, mis on siledusaste.\n")

syms x;
dfa1_dx = diff(fa1(x), x);
dfa2_dx = diff(fa2(x), x);
dfa3_dx = diff(fa3(x), x);

dfa1_dx_func = matlabFunction(dfa1_dx);
dfa2_dx_func = matlabFunction(dfa2_dx);

% Handle the constant derivative case
if isAlways(dfa3_dx == -2)
    dfa3_dx_func = @(x) -2;
else
    dfa3_dx_func = matlabFunction(dfa3_dx);
end

fprintf("Kontrollime võrdsust, et leida siledusaste.\nKontroll, kas p=1:\n")
disp(dfa1_dx_func(0.0) == dfa2_dx_func(0.0))
disp(dfa2_dx_func(1.0) == dfa3_dx_func(1.0))

fprintf("\nKontroll, kas p=2:\n")
d2fa1_dx = matlabFunction(diff(dfa1_dx_func(x), x));
d2fa2_dx = matlabFunction(diff(dfa2_dx_func(x), x));
d2fa3_dx = matlabFunction(diff(dfa3_dx_func(x), x));
disp( 2.0 == -2.0)
disp(-2.0 ==  0.0)
fprintf("Pole võrdsed, seega jääb kehtima viimane aste.\n")
fprintf("Seega on tegemist S^(2,1) splainiga.\n")


fprintf("ab)\n");
fb1 = @(x)  x;          %-1 <= x <= 0
fb2 = @(x) 1 - x;       % 0 <  x <= 1
fb3 = @(x) 2.*x - 2;    % 1 <= x <= 2

fprintf("\nKontrollime, kas esitab splaini (kui saame 0, siis ei):\n")
disp(fb1(0.0) == fb2(0.0))  % FALSE
disp(fb2(1.0) == fb3(1.0))
fprintf("Esimene võrdlus kukkus läbi, seega antud funktsioon EI ESITA splaini.\n")

%--------------------------------------------------------------------------------

fprintf("\nÜlesanne 2.\n");
x = [3.0, 4.5, 7.0, 9.0];
f = [2.5, 1.0, 2.5, 0.5];

figure(1)
plot(x, f, "o")
grid on
title("Ülesanne 2.")

S1_0 = interp1(x, f, x, "linear");
hold on
plot(x, S1_0)

val = 5.0;
S10_val = interp1(x, f, val, 'linear'); % S^(1,0) (5.0)
fprintf("S^1_0(%4.2f) = %.4f\n", val, S10_val)

%--------------------------------------------------------------------------------

fprintf("\nÜlesanne 3.\n");
f1 = @(x) 0.0;  % x <= 0
f2 = @(x) x.^2; % 0 <= x < 1
f3 = @(x) -2.*x.^2 + 6.*x - 3.0; % 1 <= x <= 2
f4 = @(x) (x-3.0).^2;   % 2 < x <= 3
f5 = @(x) 0.0;  % x >= 3

xC = [0.0, 1.0, 2.0, 3.0]; 
f = {f1, f2, f3, f4, f5};
fprintf("\nKontrollime, kas esitab splaini (kui saame 0, siis ei):\n")
for i = 1:length(xC)
    ff1 = f{i};
    ff2 = f{i+1};
    disp(ff1(xC(i)) == ff2(xC(i)))
end
syms x
df1_dx = @(x) 0.0;
df2_dx = matlabFunction(diff(f2(x), x));
df3_dx = matlabFunction(diff(f3(x), x));
df4_dx = matlabFunction(diff(f4(x), x));
df5_dx = @(x) 0.0;

fprintf("Kontrollime võrdsust, et leida siledusaste.\nKontroll, kas p=1:\n")
disp(df1_dx(0.0) == df2_dx(0.0))
disp(df2_dx(1.0) == df3_dx(1.0))
disp(df3_dx(2.0) == df4_dx(2.0))
disp(df4_dx(3.0) == df5_dx(3.0))


d2f1_dx = @(x) 0.0;
d2f2_dx = @(x) 2.0;
d2f3_dx = @(x) -4.0;
d2f4_dx = @(x) 2.0;
d2f5_dx = @(x) 0.0;

fprintf("Kontrollime võrdsust, et leida siledusaste.\nKontroll, kas p=2:\n")
disp(d2f1_dx(0.0) == d2f2_dx(0.0))
disp(d2f2_dx(1.0) == d2f3_dx(1.0))
disp(d2f3_dx(2.0) == d2f4_dx(2.0))
disp(d2f4_dx(3.0) == d2f5_dx(3.0))

%fprintf("Kolmas võrdlus kukkus läbi, seega antud funktsioon EI ESITA splaini.\n")

%--------------------------------------------------------------------------------

fprintf("\nÜlesanne 4.\n");
t = [0.0 10 20 30 40  50 60];   % min
T = [30 51 68 79 85.5 89 90];   % Celsius

tt = linspace(min(t), max(t), 500);
S3_2 = interp1(t, T, tt, "spline");
figure(2)
plot(t, T, "o")
hold on
plot(tt, S3_2)
grid on
title("Ülesanne 4.")

tval = 54; % min
S32_val = interp1(t, T, tval, 'linear'); % S^(3,2) (5.0)
fprintf("S^1_0(%4.2f) = %.4f °C\n", tval, S32_val)

%--------------------------------------------------------------------------------

fprintf("\nÜlesanne 5.\n");
Q   = [100 200 300 400 500]; % (m^3/h)
eta = [ 70  71  70  66  61]; % (%)
k = ones(1, length(Q)); % võrdsed kaalud

N = 2;  % polünoomi aste
M = N+1;
A = zeros(M,M);
B = zeros(M,1);
for i = 1:M
    B(i, 1) = sum(k.*eta.*Q.^(M - i));
    for j = 1:M
        A(i,j) = sum( k.*Q.^( 2*M - i - j ) );
    end
end
A
B
X = A^(-1) * B

% Koostame X-vektrist ruutfunkstiooni võrrandi
func = @(x) X(1).*x.^2 + X(2).*x + X(3);

figure(3)
plot(Q, eta, "o")
hold on
fplot(func, [min(Q), max(Q)])
grid on
title("Ülesanne 5.")
ylabel("Protsent %")
xlabel("m^3/h")

%--------------------------------------------------------------------------------

fprintf("\nÜlesanne 6.\n");
x = [2014 2015 2016 2017 2019 2020 2021]; % (aasta) 
y = [22.2 19.5 18.0 17.9 18.6 18.7 29.3]; % °C
k = [ 2    3    3    3    3    3    1];

% Esimese astme polünoomi koostamine
N = 1;  % polünoomi aste
M = N+1;
A1 = zeros(M,M);
B1 = zeros(M,1);
for i = 1:M
    B1(i, 1) = sum(k.*y.*x.^(M - i));
    for j = 1:M
        A1(i,j) = sum( k.*x.^( 2*M - i - j ) );
    end
end
A1
B1
X1 = A1^(-1) * B1

p1 = @(x) X1(1).*x + X1(2); 

% Teise astme polünoomi koostamine
N = 2;  % polünoomi aste
M = N+1;
A2 = zeros(M,M);
B2 = zeros(M,1);
for i = 1:M
    B2(i, 1) = sum(k.*y.*x.^(M - i));
    for j = 1:M
        A2(i,j) = sum( k.*x.^( 2*M - i - j ) );
    end
end
A2
B2
X2 = A2^(-1) * B2

p2 = @(x) X2(1).*x.^2 + X2(2).*x + X2(3);

% Ülemääratud võrrandisüsteemi abil koostatud 2. astme polünoom
N = 2;
Vtemp = vander(x);
V = Vtemp(:, (length(Vtemp) - N):length(Vtemp)) .* k';  % kaalutud Vandermonde maatriks
X_qr = linsolve(V, (k.*y)')
%[Q,R,p] = qr(V, "econ", "vector") % better. same as qr(V, 0)
%X_qr(p,:) = R\(Q \ (k .* y)')
pQR = @(x) X_qr(1).*x.^2 + X_qr(2).*x + X_qr(3);


figure(4)
plot(x, y, "o")
hold on 
fplot(p1,  [min(x), max(x)])
fplot(p2,  [min(x), max(x)])
fplot(pQR, [min(x), max(x)])
grid on
legend("Raw data", "1. astme polünoom", "2. astme polünoom", "2. astme ülemääratud võrrandisüsteem")