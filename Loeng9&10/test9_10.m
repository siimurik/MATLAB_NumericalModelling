clc
clear
fprintf("\nExercise 7\n")
x = [1.2   1.3   1.4  1.5   1.6   1.7   1.8   1.9   2];
y = [9.08 10.43 11.9 13.48 15.19 17.03 19.01 21.13 23.39];
k = [1     1     2    5     1     4     2     2     1];

figure(1)
plot(x, y, "o")
grid on

% least-squares fit (equal weights)
p = polyfit(x, y, 3)
f = @(x) p(1).*x.^3 + p(2).*x.^2 + p(3).*x + p(4);
hold on
fplot(f, [min(x), max(x)])

% least-squares fit (nonequal weights)
N = 3;      % <--- Desired polynomial degree
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

X = A\B;    % Vector for polynomial coefficients
fprintf("X = \n"); disp(X)
fM = @(x) X(1).*x.^3 + X(2).*x.^2 + X(3).*x + X(4);     % M - manual
hold on 
fplot(fM, [min(x), max(x)])

fprintf("Solution using the Vandermonde matrix and QR-decomosition\n")
% Construct the weighted Vandermonde matrix
% Starts from the last columns and moves towards the first
V = zeros(length(x), M);
for j = M:-1:1
    V(:, j) = (x.^(j - 1)) .* k;  % Scale each term by weight k
end

%Vtemp = vander(x);
%Vtemp = Vtemp(:, (length(Vtemp)-M):length(Vtemp)) .* k'
%
%Vcopy = zeros(size(Vtemp,1), size(Vtemp,2));
%for i = 1:9
%    Vcopy(i, :) = flip(Vtemp(i, :));
%end
%Vcopy
fprintf("V =\n"); disp(V)
% Solve the weighted least squares problem using QR-decomposition
% Why not use V^(-1) to solve X = V^(-1) * B? Because V is not a 
% square matrix and we cannot find a inverse for that. That is why
% we 
[Q, R] = qr(V, 0);  % not recommended but works
fprintf("Q^T =\n"); disp(Q')
fprintf("R =\n"); disp(R)
% General formula: x = inv(R) * (Q^T * (k_i ⊙ y*i)) % ⊙ -- element-wise mutiplication
X_vander = R \ (Q' * (k' .* y'));       % C(4), C(3), C(2), C(1)
x_sol=flip(X_vander);   %flipped gives us C(1), C(2), C(3), C(4)
fprintf("(flipped) X =\n"); disp(x_sol)  % C(1)*X^N + C(2)*X^(N-1) +...+ C(N)*X + C(N+1)

% Shorter way
Vtemp = vander(x);
V = Vtemp(:, (length(Vtemp) - N):length(Vtemp)) .* k'   %V - weigthed Vandermonde matrix; N - desired polynomial degree; k - weights
[Q,R,p] = qr(V, "econ", "vector") % better. same as qr(V, 0)
X_qr(p,:) = R\(Q \ (k .* y)')


%---------------------------------------------------------------------
% Solution using a special function
fprintf("Solution using the function polyfitweighted (same idea)\n")
pW = polyfitweighted(x, y, 3, k);                       % W - weighted
fprintf("pW = \n"); disp(pW)
fW = @(x) pW(1).*x.^3 + pW(2).*x.^2 + pW(3).*x + pW(4);
hold on
fplot(fW, [min(x), max(x)])
legend("Initial data", "Cubic poly (equal weights)", "Manual weighted cubic","Weighted cubic poly")

%---------------------------------------------------------------------
% A = V = S*V*D
fprintf("Solution using SVD-decomposition\n")
[L, D, R] = svd(V);
fprintf("L (left singular vectors) %dx%d\n", size(L,1), size(L,2))
disp(L)
fprintf("D (diagonal matrix of singular values) %dx%d\n", size(D,1), size(D,2)) 
disp(D)
fprintf("R (right singular vectors) %dx%d\n", size(R,1), size(R,2))
disp(R)
ky = (k .* y)';    % Weighted target vector
fprintf("Weighted target vector %dx%d\n", size(ky,1), size(ky,2))
disp(ky)

% Compute the pseudo-inverse of D aka D†
D_pseudo_inv = pinv(D);
%D_pseudo_inv = zeros(size(D'));
%for i = 1:min(size(D))
%    if D(i, i) > eps                  % Avoid division by zero
%        D_pseudo_inv(i, i) = 1 / D(i, i);
%    end
%end

% X = R ⋅ D† ⋅ L^T ⋅ ky
% Calculate the solution X
X_svd = R * D_pseudo_inv * (L' * ky);
fprintf("X_svd = \n"); disp(X_svd);

%---------------------------------------------------------------------
fprintf("Solution using Cholesky decomposition\n");
% Construct weighted matrices
W = diag(k);                % Diagonal matrix with weights
A = V' * W * V;             % Normal matrix (symmetric and positive-definite)
b = V' * W * (k.*y)';       % Weighted vector
fprintf("W (Diagonal matrix with weights) %dx%d\n", size(W,1), size(W,2))
disp(W)
fprintf("A (symmetric and positive-definite Vandermonde matrix) %dx%d\n", size(A,1), size(A,2)) 
disp(A)
fprintf("b (Weighted vector) %dx%d\n", size(b,1), size(b,2))
disp(b)

% Perform Cholesky decomposition on A
R = chol(A);                % Cholesky decomposition, A = R' * R
fprintf("Result of Cholesky decomposition A = R^T * R:\nR = \n")
disp(R)

% Solve for X in two steps using forward and back substitution
y_temp = R' \ b;            % Solve y_temp = (R^T)^-1 * b
X_chol = R \ y_temp;        % Solve X = R^-1 * y

% Display result
fprintf("Final result:\nX_chol = \n")
disp(X_chol);


%===================================================================================================
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