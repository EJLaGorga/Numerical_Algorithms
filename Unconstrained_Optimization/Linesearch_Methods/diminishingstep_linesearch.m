A = zeros(19,3);
for c = -5:2
  A(c+6,1) = (2)^c;
  [A(c+6,2), A(c+6,3)] = linesearch_decreasingstep((-2)^c);
  A(c+6, :)
end

plot(A(:,1), A(:,3))
axis([0 4 0 1.5e-5])
title('||x - min|| wrt x in decreasingstep Linesearch')
xlabel('x')
ylabel('residual of calculated min - actual min')

function [k, r] = linesearch_decreasingstep(x0)
  % parameters for linesearch with decreasing step
  x = x0';                  %initial guess of minimizer
  P = grad_f(x);            %gradient of function at guess
  k = 0;                    %counter
  
  % while the gradient is still sufficiently meaningful
  % which is equivelant to x sufficiently close to minimizer
  while  norm(P) > 1e-5  
    a = 1/(1 + k);
    P = grad_f(x);
    x = x - a*P;
    k = k + 1;
  end
    r = norm(x - [0, 0]');
end

% evaluates the gradient of (2*abs(x)^3)/3 + (X^2)/2
% at a point in R takes scalar, returns scalar
function y = grad_f(x)
  y = (2*x*abs(x)) + x;
end
