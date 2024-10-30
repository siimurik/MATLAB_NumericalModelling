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
X = A\B

A2 = vander([3 5 7]);
B2 = [-1; -1; 7];
X2 = A2\B2

hold on
fplot(@(x) X(1).*x.^2 + X(2).*x + X(3), [1, 3])
fplot(@(x) X2(1).*x.^2 + X2(2).*x + X2(3), [3, 7])

f1 = @(x) 0.0;
f2 = @(x) 2.*(x+1) + (x+1).^3;
f3 = @(x) 3.0 + 5.0.*x + 3.*x^2;
f4 = @(x) 11.0 + (x-1) + 3*(x-2).^2 + (x-1).^3;
f5 = @(x) 3.*x + 10.0;

% splain ? ... pidev üleminekupunktides
disp(f1(-1) == f2(-1))
disp(f2( 0) == f3( 0))
disp(f3( 1) == f4( 1))
disp(f4( 2) == f5( 2))
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

disp(s1(1) == s2(1))  % Seega on tegemist splainiga
% splaini järk on l = 4
% siledusaste p = ?
syms x;
s1t1 = matlabFunction(diff(s1, x, 1))
s2t1 = matlabFunction(diff(s2, x, 1))
[s1t1(1)  s2t1(1)]  % sildeusaste vähemalt  p = 1

s1t2 = matlabFunction(diff(s1t1, x, 1))
s2t2 = matlabFunction(diff(s2t1, x, 1))
[s1t2(1)  s2t2(1)]  % kukkus läbi. seega sildeusaste jääb p=1

% p = 1
% S^4,1(x)
figure(2)
fplot(s1, [0,1])
hold on
fplot(s2, [1,2])
grid on

% Ül. 4
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

% t = ?, kui s = 0.45 (miili)
interp1(s, t, 0.45, 'linear') %S^(1,0) (0.45) = 44.5200
interp1(s, t, 0.45, 'spline') %S^(3,2) (0.45) = 44.5904