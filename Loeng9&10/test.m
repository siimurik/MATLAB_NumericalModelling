clc
clear
x = [1.2   1.3   1.4  1.5   1.6   1.7   1.8   1.9   2]
y = [9.08 10.43 11.9 13.48 15.19 17.03 19.01 21.13 23.39]
k = [1     1     2    5     1     4     2     2     1]
N = 3;
Vtemp = vander(x);
V = Vtemp(:, (length(Vtemp) - N):length(Vtemp)) .* k'
[Q,R] = qr(V, 0) 
X_qr = R\(Q \ (k .* y)')

% L*D*R * X = (k .* y)'
fprintf("Lahendus SVD-teisendusega\n")
[L, D, R] = svd(V);
fprintf("L (vasakpoolsed singulaarvektorid) %dx%d\n", size(L,1), size(L,2))
disp(L)
fprintf("D (singulaarväärtuste diagonaalmaatriks) %dx%d\n", size(D,1), size(D,2)) 
disp(D)
fprintf("R (parempoolsed singulaarvektorid) %dx%d\n", size(R,1), size(R,2))
disp(R)
ky = (k .* y)';    % kaaludega parempoolne vektor
fprintf("Kaaludega parempoolne vektor %dx%d\n", size(ky,1), size(ky,2))
disp(ky)

% Pseudo-pöördmaatriks D või D†
D_pseudo_inv = pinv(D);

% Pikem viis teha sama
%D_pseudo_inv = zeros(size(D'));
%for i = 1:min(size(D))
%    if D(i, i) > eps                  % Väldi nulliga jagamist
%        D_pseudo_inv(i, i) = 1 / D(i, i);
%    end
%end

% X = R ⋅ D† ⋅ L^T ⋅ ky
% Lahendus
X_svd = R * D_pseudo_inv * (L' * ky);
fprintf("X_svd = \n"); disp(X_svd);

fprintf("Lahendus Cholesky teisenduse abil\n");
% Koosta kaalutud maatriksid
W = diag(k);                % Diagonaalmaatriks kaaludega
A = V' * W * V;             % Normaalmaatriks (sümmeetriline ja positiivne definitsiooniga)
b = V' * W * (k.*y)';       % Kaalutud vektor
fprintf("W (Diagonaalmaatriks kaaludega) %dx%d\n", size(W,1), size(W,2))
disp(W)
fprintf("A (sümmeetriline ja positiivse Vandermonde maatriks) %dx%d\n", size(A,1), size(A,2)) 
disp(A)
fprintf("b (Kaalutud vektor) %dx%d\n", size(b,1), size(b,2))
disp(b)

% Cholesky teisendus maatriksile A
R = chol(A);                % Cholesky teisendus, A = R' * R
fprintf("Cholesky teisenduse tulemus A = R^T * R:\nR = \n")
disp(R)

% Lahendus kahes etapis: edasi- ja tagasisubstitutsioon
y_temp = R' \ b;            % Lahenda y_temp = (R^T)^-1 * b
X_chol = R \ y_temp;        % Lahenda X = R^-1 * y

% Kuva tulemus
fprintf("Lõpptulemus:\nX_chol = \n")
disp(X_chol);