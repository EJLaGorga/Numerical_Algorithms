rand('seed', 314);
x = linspace(0,3,30);
y = 2*x.^2-3*x+1+5*rand(size(x));


X = zeros(30, 2);
X(:,1) = 1;
X(:,2) = x(:);
X(:,3) = x(:).^2;

% Normal Equations
p = inv(X'*X)*X'*y';

% QR Factorization
[Q,R] = qr(X);
p = R\Q'*y';

% SVD Decomposition
[U,E,V] = svd(X, 'econ');
p = V*inv(E)*U'*y';

scatter(x, y);
hold on
fplot(@(x) p(1) + p(2)*x + p(3)*x.^2);