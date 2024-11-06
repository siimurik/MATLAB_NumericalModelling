clc
clear

t = [5, 8, 10, 12, 15, 19, 22, 25]
v = [106.8, 177.2, 232.6, 279.2, 301.4, 372.3, 401.3, 442.5]
figure(1)
plot(t, v, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
grid on
a = 2.0;
b = 4.0;
c = 5.0;
d = 8.0;
e = 4.0;
f = 1.0;
k = [1, (1+a), (1+b), (1+a), (1+c), (1+d), (1+e), (1+f)]

fprintf("a) Lineaarfunktsiooniga vähimruutude mõttes võrdsete kaaludega")
% lineaarfunktsiooniga vähimruutude mõttes
p1 = polyfit(t, v, 1)    % lineaarfuktsioon
f1 = @(x) p1(1).*x + p1(2);

% Ruutfunktiooniga
p2 = polyfit(t, v, 2)    % ruutfunktsioon
f2 = @(x) p2(1).*x.^2 + p2(2).*x + p2(3);
hold on
fplot(f1, 'm')
fplot(f2, 'b')
title("Võrdsete kaaludega funktsioonid")
legend("Raw data", "f1(x) = 16.1010*x + 55.6978", "f2(x) = -0.4198*x^2 + 28.8003*x - 22.0144", "Location", "southeast")


fprintf("b) Lineaarfunktsiooniga vähimruutude mõttes etteantud kaaludega")
N = 1;  % polünoomi aste
M = N+1;
A_w1 = zeros(M,M);
B_w1 = zeros(M,1);
for i = 1:M
    B_w1(i, 1) = sum(k.*v.*t.^(M - i));
    for j = 1:M
        A_w1(i,j) = sum( k.*t.^( 2*M - i - j ) );
    end
end
A_w1
B_w1
X_w1 = A_w1^(-1) * B_w1
fw1 = @(x) X_w1(1).*x + X_w1(2);

N = 2;  % polünoomi aste
M = N+1;
A_w2 = zeros(M,M);
B_w2 = zeros(M,1);
for i = 1:M
    B_w2(i, 1) = sum(k.*v.*t.^(M - i));
    for j = 1:M
        A_w2(i,j) = sum( k.*t.^( 2*M - i - j ) );
    end
end
A_w2
B_w2
X_w2 = A_w2^(-1) * B_w2
fw2 = @(x) X_w2(1).*x.^2 + X_w2(2).*x + X_w2(3);

figure(2)
plot(t, v, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
grid on
hold on
fplot(fw1, 'm')
fplot(fw2, 'b')
legend("Raw data", "fw1(x) = 15.5011*x + 70.0295", "fw2(x) = -0.3632*x^2 + 26.6967*x - 5.8654", "Location", "southeast")
title("Erikaaludega funktsioonid")

fprintf("c) Splainidega lähendamine\n")
figure(3)
plot(t,v,'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
xlabel('s')
ylabel('t')
grid on

S1_0 = interp1(t, v, t, "linear");
hold on
plot(t, S1_0, 'b')

tt = [min(t): 0.001: max(t)];
S3_2 = interp1(t, v, tt, "spline");
hold on
plot(tt, S3_2, 'm')
hold off
title("Splainidega graafik")
legend("Raw data", "S^{1,0}(x) splain", "S^{3,2}(x) splain", "Location", "southeast")

fprintf("d) raketi kiirus ajahetkel 17 sekundit?\n")

fprintf("Võrdsete kaaludega lineaarfunktsioon vähimruutude mõttes: f1(17.0) = %.4f s.\n", f1(17.0))
fprintf("Võrdsete kaaludega ruutfunktsioon vähimruutude mõttes: f2(17.0) = %.4f s.\n",    f2(17.0))
fprintf("Erinevate kaaludega lineaarfunktsioon vähimruutude mõttes: fw1(17.0) = %.4f s.\n",fw1(17.0))
fprintf("Erinevate kaaludega ruutfunktsioon vähimruutude mõttes: fw2(17.0) = %.4f s.\n",   fw2(17.0))
S10_val = interp1(t, v, 17.0, 'linear'); 
S32_val = interp1(t, v, 17.0, 'spline'); 
fprintf("Lineaarsplain: S^1_0(%4.2f) = %.4f s.\n", 17.0, S10_val)
fprintf("Kuupsplain:    S^3_2(%4.2f) = %.4f s.\n", 17.0, S32_val)