function [A] = gen_lmu(s,m,n)
%GENLMU Summary of this function goes here
%   Detailed explanation goes here 
if nargin < 3
    U = orth(rand(m));
    e = 10*rand(1,m);
    e(1,1) = s;
    e(1,m) = 1/s;
    E = diag(e);
    V = orth(rand(m));
    A = U*E*V';
else

if m < n
    [m, n] = deal(n,m);
end
U = orth(rand(m));
e = 10*rand(1,n);
e(1,1) = s;
e(1,n) = 1;
E = diag(e, m-n);
E = E(:, (m-n+1):end);
V = orth(rand(n));
A = U*E*V';
end
end