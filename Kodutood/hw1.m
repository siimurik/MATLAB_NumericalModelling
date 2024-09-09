clc
clear
% Iseseisev töö nr. 1

disp("Ülesanne 1.")
fprintf('a) %.4e\n', 2^(4^3)        );
fprintf('b) %.4f\n', (2^4)^3        );
fprintf('c) %.4f\n', log(cos(10^0)) );
fprintf('d) %.4f\n\n', log10(-pi)   );

disp("Ülesanne 2.")
fprintf('a) %.4f\n', 5*sqrt(3)/8 + 15*(1/13 - 2/(5*4^3) - 1/(8*2^5)) );
fprintf('b) %.4f\n\n', sqrt(2 + sqrt(2 + sqrt(2 + sqrt(2)))) );

disp("Ülesanne 3.")
tic;
(-64.)^(1/3)    % NOTE: What????
toc

tic;
-64.^(1/3)      % Normaalne viis. Kõige kiirem.
toc

tic;
nthroot(-64, 3)
toc
fprintf("\n");

disp("Ülesanne 4.")
r = 3;
h = 13;
V = pi*(r^2)*h
fprintf("\n");

disp("Ülesanne 5.")
x4  = 4.0
mn_ = 2.0
pi  = 3
%_x7 = 7.0
y2z = 2.0
%3y = 3.0
fprintf("\n");

disp("Ülesanne 6.")
[1:6]
[5:9]
[1:3:10]
[3:2:14]
fprintf("\n");

disp("Ülesanne 7.")
disp("y = linspace(x1,x2) tagastab rea vektori, mis koosneb 100 ühtlaselt jaotatud punktist vahemikus x1 kuni x2.")
disp("y = linspace(x1,x2,n) genereerib n punkti. Punktide vahe on (x2-x1)/(n-1).")
fprintf("\n");

disp("Ülesanne 8.")
linspace(1,4,2)
linspace(1,4,6)
fprintf("\n");

disp("Ülesanne 9.")
disp("a ^ 3  - astendamine")
disp("a .^ 3 - komponentide kaupa astendamine")
fprintf("\n");

disp("Ülesanne 10.")
sum = 0.0;
for i = 1:51    % 2n-1=101 --> 2n=102 --> n=51
    sum = sum + 1/(2*i-1);
end
sum
fprintf("\n");

disp("Ülesanne 11.")
alpha = [0., 30., 45., 60., 90.];
for j = 1:length(alpha)
    disp( (sin(alpha(j)))^2 + (cos(alpha(j)))^2 == 1 )
end
fprintf("\n");

disp("Ülesanne 12.")
A = [ 4. -1. 3. 10.
      0.  7. 9. 12.
     21. -6. 2.  8.];

B = [15. -12.  0.  4.
      1.   6. -3.  1.
     -5.   5.  14.  20.]; 

C = [-4.   7.  3. 11.
     31. -18. -2.  4.
     70.   1.  2. 22.];

(5*A - 8*B) * transpose(C)
fprintf("\n");

disp("Ülesanne 13.")
[0, 1; 2, 3]^2
[0, 1; 2, 3].^2
fprintf("\n");

disp("Ülesanne 14.")
clear A B 
A = [1. -1.  3. 0.
     2.  0. -1. 2.];

B = [4. 2. -1. 2.
     3. 1.  5. 0.];
    
A+B
A.*B
A/B
A./B
A.^B
fprintf("\n");

disp("Ülesanne 15.")
clear A B
A = [1  5; 12 7; 4 3]
B = [-9 2 9; 15 6 11]

%C = A*B;
%fprintf('a) C = A*B %.4f\n', C);
%fprintf('b) C^T, C^-1, det(C) %.4f\n', C);

disp("a) C = A*B")
C = A*B

disp("b) C^T")
transpose(C)

disp("C^-1")
inv(C)    % NOTE: Pöördmaatriksit vist ei leidu.

disp("det(C)")
det(C)

fprintf("\nc) a_22 - 6b_13\n");
A(2,2) - 6*B(1,3)
fprintf("\n");


disp("Ülesanne 16.")
a = [4, 3 , -5, 6];
b = [8; -2; 1; 9];
a*b
fprintf("\n");

disp("Ülesanne 17.")
clear A
A = [1 2 3 4 5 6 7; 9 7 5 3 1 -1 -3; 4 8 16 32 64 128 256]

%Alternatiivselt
A = zeros(3, 7);
for n = 1:7
    A(1, n) = n;
    A(2, n) = 9-2*(n-1);
    A(3, n) = 2^(2+(n-1));
end
A
fprintf("\n");

disp("Ülesanne 18.")
clear A
A = [ 5.  6. -2. 10.  4.  7.
      8.  1.  3.  0. -9. 12.
     -6.  5.  7. 24. 41.  5.
      3. -3. 15. -6.  2. -7.
      0. 41.  9. -1.  3.  8.];

A(4, 3)
A(4, :)
A(:, 2)
fprintf("\n");

disp("Ülesanne 19.")
a = rand(1, 10);
b = rand(1, 10);

disp('a)')
disp(a.^3 .* sqrt(b))

disp('b)')
disp( nthroot(a, 5) ./ (b - 3*log(a)) )
fprintf("\n");

disp("Ülesanne 20.")
disp("I = eye(n) tagastab n x n ühikmaatriksi, mille peadiagonaalil on ühed ja mujal nullid.")
disp("Suurus n on maatriksi ridade ja veergude arv. Samuti ka ühtede arv peadiagonaalil.")
fprintf("\n");

disp("Ülesanne 21.")
y = @(x) 7.*x.^4 + 3.*x.^3 + 4.*x.^2 - 2.*x + 8   % warning-free version
%y = @(x) 7*x^4 + 3*x^3 + 4*x^2 - 2*x + 8         % "improper vectorization"
figure(3)
fplot(y, [-5, 5], 'green')

disp("Ülesanne 22.")
x = [4 2 7]
y = [-5 4 3]
figure(1)
plot(x, y)
xlabel('x-axis')
ylabel('y-axis')
title('Plot of x vs y')

disp("Ülesanne 23.")
clear pi
pi = 4.0*atan(1.0);
x = 0:pi/100:2*pi;
y = @(x) x .* exp(sin(x)); 
figure(2)
fplot(y, [-3, 2], 'red')
