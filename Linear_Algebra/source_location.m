randn('seed', 317);
% each c||olumn of A represents lat lon of sensor
A = randn (2,5);
% location of source used for data generation
x = randn (2,1);
% noisy measurements of distances between source and sensors
d = sqrt(sum((A-x*ones(1,5)).^2))+0.05*randn(1,5);

% I finally learned how to use anonomous funtions in matlab!
residual = @(x) (sum(([x x x x x] - A).^2, 1) - d.^2)';
jacobian = @(x) 2*([x x x x x] - A)';
evaluate = @(x) norm(residual(x))^2;

p = .1;
c = .5;
k = 0;
xk = [1;1];
XK = [];
Pk = 1;

% Gauss-Newton Method for nonlinear least squares
while norm(Pk) > 1e-6
    Jk = jacobian(xk);          % (5x2)
    rk = residual(xk);          % (5x1)
    Pk = - inv(Jk'*Jk)*Jk'*rk;  % (2x1)
    
    fk = evaluate (xk);
    a = 1;
    
    while evaluate(xk+a*Pk) > fk+c*a*((Jk'*rk)'*Pk)
        a = p*a;
    end
    
    XK = [XK, norm(xk + a*Pk - x)/norm(xk - x)];
    xk = xk + a*Pk
    k = k + 1;
end

XK'