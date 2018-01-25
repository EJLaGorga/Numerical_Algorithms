function [KL] = approxKLdiv(Mu, Sigma, MuR, SigmaR)
% Approximate KL divergence of guassian mixture
% by using equiweighted sums and closed form 
% KL divergence for induvidual gaussians

assert(size(Mu,1) == size(Sigma,3) && size(Mu,1) == size(MuR,1));
k = size(Mu,1);
n = size(Mu,2);

mu0 = sum(Mu,1);
sigma0 = zeros(n,n);
for i = 1:k
    sigma0 = Sigma(:,:,i) + (Mu(i,:) - mu0)*(Mu(i,:) - mu0)';
end

mu1 = sum(MuR,1);
sigma1 = zeros(n,n);
for i = 1:k
    sigma1 = SigmaR(:,:,i) + (MuR(i,:) - mu1)*(MuR(i,:) - mu1)';
end

tmp = inv(sigma1)*sigma0;
KL = 0.5*(trace(tmp)+(mu1-mu0)*inv(sigma1)*(mu1-mu0)'-n-log(det(tmp)));