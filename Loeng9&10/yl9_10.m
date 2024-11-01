clc
clear

x = [1, 2,  3,  5, 7];
y = [5, 1, -1, -1, 7];
figure(1)
plot(x, y, "o")
grid on
axis([-2 9 -3 9])

A = vander([1 2 3]);
B = [5; 1; -1];
X = A\B;
fprintf("Võrrandisüsteemi lahendid:\n X ="); disp(X)

A2 = vander([3 5 7]);
B2 = [-1; -1; 7];
X2 = A2\B2;
fprintf("Võrrandisüsteemi lahendid:\n X2 ="); disp(X2)

hold on
fplot(@(x) X(1).*x.^2 + X(2).*x + X(3), [1, 3])
fplot(@(x) X2(1).*x.^2 + X2(2).*x + X2(3), [3, 7])

f1 = @(x) 0.0;
f2 = @(x) 2.*(x+1) + (x+1).^3;
f3 = @(x) 3.0 + 5.0.*x + 3.*x^2;
f4 = @(x) 11.0 + (x-1) + 3*(x-2).^2 + (x-1).^3;
f5 = @(x) 3.*x + 10.0;

% splain ? ... pidev üleminekupunktides
fprintf("Kontrollime, kas esitab splaini (kui saame 0, siis ei):\n")
fprintf("%d\n", f1(-1) == f2(-1))
fprintf("%d\n", f2( 0) == f3( 0))
fprintf("%d\n", f3( 1) == f4( 1))
fprintf("%d\n", f4( 2) == f5( 2))
% seega f(x) ei esita splaini

%{
f = {f1, f2, f3, f4, f5};
xPoints = [-1, 0, 1, 2];
for i=1:length(f)-1
    fn  = f(i);
    fn1 = f(i+1);
    fn(xPoints(i)) == fn1(xPoints(i))
end
%}

% Ül 3
clear
s1 = @(x) x.^4 + - 2.*x.^3 - 2.*x.^2 - 4.*x - 6;
s2 = @(x) x.^4 + - 2.*x.^3 - 4.*x.^2 - 8;

fprintf("\nKontrollime, kas esitab splaini (kui saame 0, siis ei):\n")
disp(s1(1) == s2(1))  % Seega on tegemist splainiga
% splaini järk on l = 4
% siledusaste p = ?
syms x;
s1t1 = matlabFunction(diff(s1, x, 1));
s2t1 = matlabFunction(diff(s2, x, 1));
fprintf("Kontrollime võrdsust, et leida siledusaste.\nKontroll, kas p=1:\n")
disp([s1t1(1)  s2t1(1)])  % sildeusaste vähemalt  p = 1

s1t2 = matlabFunction(diff(s1t1, x, 1));
s2t2 = matlabFunction(diff(s2t1, x, 1));
fprintf("Kontrollime võrdsust, et leida siledusaste.\nKontroll, kas p=2:\n")
disp([s1t2(1)  s2t2(1)])  % kukkus läbi. seega sildeusaste jääb p=1
fprintf("Pole võrdsed, seega jääb kehtima viimane aste.\n")

% p = 1
% S^4,1(x)
figure(2)
fplot(s1, [0,1])
hold on
fplot(s2, [1,2])
grid on

% Ül. 4
fprintf("\nÜlesanne 4.\n")
clear
s = [0, 0.25,  0.5, 0.75, 1.0,  1.25];
t = [0,   25, 49.4,   73, 96.4, 119.4];
figure(3)
plot(s,t,'o')
xlabel('s')
ylabel('t')
grid on

S1_0 = interp1(s, t, s, "linear");
hold on
plot(s, S1_0, 'b')

ss = [0: 0.001: 1.25];
S3_2 = interp1(s, t, ss, "spline");
hold on
plot(ss, S3_2, 'r')
hold off

% t = ?, kui s = 0.45 (miili)
S10_val = interp1(s, t, 0.45, 'linear'); %S^(1,0) (0.45) = 44.5200
S32_val = interp1(s, t, 0.45, 'spline'); %S^(3,2) (0.45) = 44.5904
fprintf("S^1_0(%4.2f) = %.4f\n", 0.45, S10_val)
fprintf("S^3_2(%4.2f) = %.4f\n", 0.45, S32_val)

fprintf("\nÜlesanne 5\n")
clear
x = [1 2  5  8 11 12]
y = [5 8 14 26 50 98]
figure(4)
plot(x, y, "o")
grid on
title("Ylesanne 5")


% PIKEMALT: Meetod 1
% Võrdsed kaalud, ehk kappa = 1 
% Võib ka olla midagi muud sest nagunii taandub lõpus kõik välja
k = ones(1, length(x))

% lineaarfunktsiooniga vähimruutude mõttes
% Phi(x) = c1*x + c0
A = [   sum(x.^2),      sum(x)      % sum = 0.0; for i=1:length(x)
        sum(x)          sum(k)]     %       sum = sum + x(i)*x(i);
                                    % end; sum
% Üldine meetod
N = 1; M = N+1;
A_copy = zeros(M,M);
B_copy = zeros(M,1);
for i = 1:M
    B_copy(i, 1) = sum(k.*y.*x.^(M - i));
    for j = 1:M
        A_copy(i,j) = sum( k.*x.^( 2*M - i - j ) );
    end
end
A_copy

B = [   sum(y.*x)
        sum(y)]
B_copy

X = A^(-1) * B

f1 = @(x) X(1).*x + X(2);
% Phi(x) ~ 6.8104*x - 10.7678
hold on
fplot(f1)

% Meetod 2: lühemalt
p1 = polyfit(x, y, 1)    % leiab konstantid ainult etteantud x ja y vektorite põhjal
xx = [3, 9];
polyval(p1, xx)

% 2. ruutfunktsiooniga vähimruutude mõttes
p2 = polyfit(x, y, 2)
fancyFunc = polyfitweighted(x,y,2,k)

f2 = @(x) p2(1).*x.^2 + p2(2).*x + p2(3);
hold on
fplot(f2)
polyval(p2, xx)

fprintf("\nÜlesanne 7\n")
clear
x = [1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2];
y = [9.08 10.43 11.9 13.48 15.19 17.03 19.01 21.13 23.39];
k = [1 1 2 5 1 4 2 2 1];

% kuupsplain
xx = [min(x): 0.001: max(x)];
S3_2 = interp1(x, y, xx, "spline");
figure(5)
plot(x, y, "o")
hold on
plot(xx, S3_2, 'r')
grid on
title("Ylesanne 7")

S32_X1 = interp1(x, y, 1.12, "spline");
S32_X2 = interp1(x, y, 1.35, "spline");

x_val = [1.12, 1.35];
for i = 1:length(x_val)
    fprintf("S3_2(%3.2f) = %f\n", x_val(i), interp1(x, y, x_val(i), "spline"))
end

% vähimruutude meetod (võrdsed kaalud)
p = polyfit(x, y, 3)
f = @(x) p(1).*x.^3 + p(2).*x.^2 + p(3).*x + p(4);
hold on
fplot(f, [min(x), max(x)])

% vähimruutude meetod (mittevõrdsed kaalud)
N = 3;      % Soovitud polünoomi aste
M = N+1;    % "4"
A = zeros(M, M);
X = zeros(M, M);
B = zeros(M, 1);

for i = 1:length(X)
    B(i, 1) = sum(k.*y.*x.^(M - i));
    for j = 1:length(X)
        A(i, j) = sum( k.*x.^( 2*M - i - j ) );   % 8 = 2*length(X)???
    end
end
fprintf("A = \n"); disp(A)
fprintf("B = \n"); disp(B)

X = A\B;
fprintf("X = \n"); disp(X)
fM = @(x) X(1).*x.^3 + X(2).*x.^2 + X(3).*x + X(4);     % M - manual (käsitsi koostatud)
hold on 
fplot(fM, [min(x), max(x)])

pW = polyfitweighted(x, y, 3, k);                       % W - weighted
fprintf("pW = \n"); disp(pW)
fW = @(x) pW(1).*x.^3 + pW(2).*x.^2 + pW(3).*x + pW(4);
hold on
fplot(fW, [min(x), max(x)])
legend("Initial data", "S^3_2 spline", "Cubic poly (equal weights)", "Manual weighted cubic","Weighted cubic poly")

function p = polyfitweighted(x,y,n,w)
    % polyfitweighted.m 
    % -----------------
    %
    % Find a least-squares fit of 1D data y(x) with an nth order 
    % polynomial, weighted by w(x).
    %
    % By S.S. Rogers (2006), based on polyfit.m by The MathWorks, Inc. - see doc
    % polyfit for more details.
    %
    % Usage
    % -----
    %
    % P = polyfitweighted(X,Y,N,W) finds the coefficients of a polynomial 
    % P(X) of degree N that fits the data Y best in a least-squares sense. P 
    % is a row vector of length N+1 containing the polynomial coefficients in
    % descending powers, P(1)*X^N + P(2)*X^(N-1) +...+ P(N)*X + P(N+1). W is
    % a vector of weights. 
    %
    % Vectors X,Y,W must be the same length.
    %
    % Class support for inputs X,Y,W:
    %    float: double, single
    %
    % The regression problem is formulated in matrix format as:
    %
    %    yw = V*p    or
    %
    %          3    2
    %    yw = [x w  x w  xw  w] [p3
    %                            p2
    %                            p1
    %                            p0]
    %
    % where the vector p contains the coefficients to be found.  For a
    % 7th order polynomial, matrix V would be:
    %
    % V = [w.*x.^7 w.*x.^6 w.*x.^5 w.*x.^4 w.*x.^3 w.*x.^2 w.*x w];
    if ~isequal(size(x),size(y),size(w))
        error('X and Y vectors must be the same size.')
    end
    x = x(:);
    y = y(:);
    w = w(:);
    % Construct weighted Vandermonde matrix.
    %V(:,n+1) = ones(length(x),1,class(x));
    V(:,n+1) = w;
    for j = n:-1:1
       V(:,j) = x.*V(:,j+1);
    end
    % Solve least squares problem.
    [Q,R] = qr(V,0);
    ws = warning('off','all'); 
    p = R\(Q'*(w.*y));    % Same as p = V\(w.*y);
    warning(ws);
    if size(R,2) > size(R,1)
       warning('polyfitweighted:PolyNotUnique', ...
           'Polynomial is not unique; degree >= number of data points.')
    elseif condest(R) > 1.0e10
        if nargout > 2
            warning('polyfitweighted:RepeatedPoints', ...
                'Polynomial is badly conditioned. Remove repeated data points.')
        else
            warning('polyfitweighted:RepeatedPointsOrRescale', ...
                ['Polynomial is badly conditioned. Remove repeated data points\n' ...
                '         or try centering and scaling as described in HELP POLYFIT.'])
        end
    end
    p = p.';          % Polynomial coefficients are row vectors by convention.
end