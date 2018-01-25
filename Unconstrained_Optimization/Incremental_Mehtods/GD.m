function [nXk] = GD(xk)
%IGD Summary of this function goes here
%   Detailed explanation goes here
m = 10;
n = 2;
eta = 1;
[A, b] = generate_Ab(m, n, eta);
gradient = @(x) gradient_f(A, b, x);

nXk = [];
x = [0.043050527109744; 0.485640581896102];
y = 1; 
L = 100;
a = 2/(y+L);

for k = 0:499
    nXk = [nXk; log(norm(xk - x))];
    Pk = gradient(xk);
    xk = xk - a*Pk;
end

end