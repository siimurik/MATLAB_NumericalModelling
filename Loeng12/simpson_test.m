clc
clear
a = pi/4;
b = pi/2;
n = 10;
if mod(n,2) == 0    % kui jääk on 0, siis n on paarisarv
    h = (b-a)/n;
else                % paaritu arvu korral muuda osalõikude arv ühe võrra suuremaks 
    n = n+1;
    h = (b-a)/n;
end
x = linspace(a, b, n+1); % n+1 sõlme n osalõigu jaoks
f = @(x) sin(x)./x;

I_simp = (h/3) * (f(x(1)) + 2*sum(f(x(3:2:end-2))) + 4*sum(f(x(2:2:end-1))) + f(x(end)));
I_exact = integral(f, a, b);
err_simp = I_exact - I_simp;

disp(['Simpsoni reegli lähendus: ', num2str(I_simp)])
disp(['Viga: ', num2str(err_simp)])
disp(['MATLAB''i integraalfunktsiooni tulemus: ', num2str(I_exact)])