function [ pdf ] = gaussianND(X, Mu, Sigma)
% Calculate multivariate Gaussian
% X     - Vector of data points
% Mu    - Mean Vector
% Sigma - Covariance matrix

n = size(X, 2);
meanDiff = bsxfun(@minus, X, Mu);
pdf = 1 / sqrt((2*pi)^n * det(Sigma)) * ...
      exp(-1/2 * sum((meanDiff * inv(Sigma) .* meanDiff), 2));
end
