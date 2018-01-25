function [X] = gen2DGaussianData(Mu, Sigma, Ns)
% Generate n data points from guassian distribution
% with the parameters mu, sigma. Plot data and
% contour of probability density function.
% Return mixed data set containing all data points.
% Mu    - Vector of mean for each distribution
% Sigma - Vector of covarience matricies
% Ns    - Vector of number of data to be generate

assert(size(Mu,1) == size(Sigma,3) && size(Mu,1) == size(Ns,1))
%assert(size(Mu,2) == 2)

X = [];
figure(1);
% Generate and plot sample data from the distributions
for i = 1:size(Mu,1)
    Xs = randn(Ns(i),size(Mu,2))*chol(Sigma(:,:,i));
    Xs = Xs + repmat(Mu(i,:), Ns(i), 1);
    X = [X; Xs];
    plot(Xs(:, 1), Xs(:, 2), 'o');
    hold on;
end


if size(Mu,2) == 2
% Plot the contour lines of the distributions
gridSize = 100;
Mx = max(X(:,1))+1;
mx = min(X(:,1))-1;
My = max(X(:,2))+1;
my = min(X(:,2))-1;

u = linspace(min([mx,my]), max([Mx,My]), gridSize);
[A, B] = meshgrid(u, u);
gridX = [A(:), B(:)];

for i = 1:size(Mu,1)
    pd = gaussianND(gridX, Mu(i,:), Sigma(:,:,i));
    PD = reshape(pd, gridSize, gridSize);
    contour(u, u, PD);
end

axis([mx,Mx,my,My])
title('Generated Data and Actual Distributions');
end
end
