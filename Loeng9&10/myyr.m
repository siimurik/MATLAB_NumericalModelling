clc
clear

t = [5, 8, 10, 12, 15, 19, 22, 25]
v = [106.8, 177.2, 232.6, 279.2, 301.4, 372.3, 401.3, 442.5]
k = [1 3 3 3 2 4 8 6];

% Lineaarfunktsiooni korral
A1 = [  sum(k.*t.^2),    sum(k.*t.^1);      
        sum(k.*t.^1),    sum(k.*t.^0)]

B1 = [  sum(k.*v.*t.^1)
        sum(k.*v.*t.^0)]
X1 = A1^(-1) * B1

% Ruutfunktsiooni korral
A2 = [  sum(k.*t.^4), sum(k.*t.^3),    sum(k.*t.^2);      
        sum(k.*t.^3), sum(k.*t.^2),    sum(k.*t.^1);
        sum(k.*t.^2), sum(k.*t.^1),    sum(k.*t.^0)]

B2 = [  sum(k.*v.*t.^2)
        sum(k.*v.*t.^1)
        sum(k.*v.*t.^0)]
X2 = A2^(-1) * B2

% Meetod 2:
N = 2;  % Lineaarfunkt - "1"; ruutfunkt - "2"
M = N+1;
A = zeros(M,M);
B = zeros(M,1);
for i = 1:M
    B(i, 1) = sum(k.*v.*t.^(M - i));
    for j = 1:M
        A(i,j) = sum( k.*t.^( 2*M - i - j ) );
    end
end
A
B
X = A \ B

% Meetod 3: Most advanced and accurate
N = 2;  % Lineaarfunkt - "1"; ruutfunkt - "2"
Vtemp = vander(t);
V = Vtemp(:, (length(Vtemp) - N):length(Vtemp)) .* k' 
[Q,R,p] = qr(V, "econ", "vector") % better. same as qr(V, 0)
X_qr(p,:) = R\(Q \ (k .* v)')
