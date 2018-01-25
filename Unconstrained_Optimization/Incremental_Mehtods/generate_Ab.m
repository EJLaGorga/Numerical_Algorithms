function [ A, b ] = generate_Ab(m, n, eta)
% generate 
%
% @ Xiaozhe Hu

% initialize
A = cell(m,1);
b = cell(m,1);

% fix seed
rng(1); 

for i = 1:m
    
    % generate A
    diag_A = [10.^randi([0 eta],n/2,1); 10.^randi([-eta 0],n/2,1)];
    A{i} = spdiags(diag_A, 0, n, n);
    
    % generate b
    b{i} = rand(n,1);
    
end

end