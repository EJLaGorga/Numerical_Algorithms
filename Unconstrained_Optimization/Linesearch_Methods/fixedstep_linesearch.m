A = zeros(19,3);
for c = 1:19
  A(c,1) = c/20;
  [A(c,2), A(c,3)] = linesearch_fixedstep(c/20);
end

%figure
%subplot(2,1,1)
%plot(A(:,1), A(:,2))
%title('Convergence wrt \alpha in fixedstep Linesearch')
%xlabel('\alpha')
%ylabel('Iterations before ||-\Deltaf(x)|| <= 10^{-5}')

%subplot(2,1,2)
plot(A(:,1), A(:,3))
axis([0 1 1e-6 10e-6])
title('||x - min|| wrt \alpha in fixedstep Linesearch')
xlabel('\alpha')
ylabel('residual of calculated min - actual min')

function [k, r] = linesearch_fixedstep(q)
  % parameters for linesearch with fixedstep
  x = [0, 0]';              %initial guess of minimizer
  P = grad_f(x);            %gradient of function at guess
  a = q;                    %random scalar in (0,1)
  k = 0;                    %counter
  
  % while the gradient is still sufficiently meaningful
  % which is equivelant to x sufficiently close to minimizer
  while  norm(P) > 1e-5  
    P = grad_f(x);
    x = x - a*P;
    k = k + 1;
  end
  r = norm(x - [2, 1]');
end

% evaluates the gradient of (1/2)x'Ax - b'x
% at a point in R^2 takes 2d vector, returns 2d vector
function y = grad_f(x)
  y(1) = 2*x(1) - x(2) - 3;
  y(2) = 2*x(2) - x(1);
  y = y';
end
