clc
clear

% Ülesanne 1
f = @(x) x.^3 - 8.*x + 2;
figure(1)
fplot(f)
yline(0,"k-") % y = 0 ehk x-telg
xline(0,"k-") % x = 0 ehk y-telg
grid on
hold off
% määrame võrrandi alglähendid jooniselt
% võtame alglähenditeks: -3, 0.2, 2.5
f(-3)
lahend1 = fzero(f, -3)
kontroll1 = f(lahend1)

lahend2 = fzero(f, 0.2)
kontroll2 = f(lahend2)

lahend3 = fzero(f, 2.5)
kontroll3 = f(lahend3)

% kui mitu lahendit on piirkonnas, läheb fzero() segadusse
%fzero(f, [-4, 1])   

% Ülesanne 2
y = @(x) log(x) + x;
figure(2)
fplot(y)
yline(0,"k-") % y = 0 ehk x-telg
xline(0,"k-") % x = 0 ehk y-telg
grid on
hold off

lahend = fzero(y, [0.1, 1.0])

% Ülesanne 3
g = @(x) x.^3 - 6.*x.^2 + 3.*atan(x) + 3;
figure(3)
fplot(g)
yline(0,"k-") % y = 0 ehk x-telg
xline(0,"k-") % x = 0 ehk y-telg
grid on
hold off

lahend1 = fzero(g, -0.5)
kontroll1 = g(lahend1)

lahend2 = fzero(g, 1)
kontroll2 = g(lahend2)

lahend3 = fzero(g, 5.5)
kontroll3 = g(lahend3)

% Ülesanne 4
% võrrandi lahendamine hariliku iteratsioonimeetodiga (HIM)
f = @(x) x.^3 - 8.*x + 2;
figure(4)
fplot(f)
yline(0,"k-") % y = 0 ehk x-telg
xline(0,"k-") % x = 0 ehk y-telg
grid on
hold off
% 1) alglähendid määrame graafiliselt
% alglähendid: -3, 0.2, ja 2.5 .

% 2) avaldame võrrandist f(x) = 0 suuruse x, saame x=g(x)
g = @(x) 1./8 .* (x.^3 + 2);
dgdx = @(x) (3.*x.^2)./8;
c = [-3, 0.2, 2.5];
abs(dgdx(c)) < 1    % 0 = ei sobi; 1 = sobib; 0 = ei sobi

% syms x
% g = nthroot(8*x - 2, 3);
% dgdx = diff(g)
x = 0.2;
count = 0;
while abs(f(x)) > 1.0E-6
    x = g(x);
    count = count + 1;
end
lahend2 = x
sammud2 = count

% Proovime uuesti
g = @(x) nthroot(8.*x - 2, 3);
dgdx = @(x) 8./(3.*(8.*x - 2).^(2/3));
abs(dgdx(c)) < 1    % 1 = sobib; 0 = ei sobi; 1 = sobib

x = c(1)
count = 0;
while abs(f(x)) > 1.0E-6
    x = g(x);
    count = count + 1;
end
lahend1 = x
sammud1 = count

x = c(3)
count = 0;
while abs(f(x)) > 1.0E-6
    x = g(x);
    count = count + 1;
end
lahend3 = x
sammud3 = count