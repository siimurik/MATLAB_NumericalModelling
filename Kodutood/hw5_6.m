clc
clear


% Ylesanne 2
A = [ 8.0  7.0 -2.0
     -6.0  4.0  3.0
      4.0 -1.0  5.0];

norm(A, 1)
norm(A, inf)

% Ylesanne 3
F = @(z) [ z(1)^3 + z(2) - 1, z(2)^3 - z(1) + 1];
x0 = [0, 0];
z_lahend = fsolve(F, x0)
