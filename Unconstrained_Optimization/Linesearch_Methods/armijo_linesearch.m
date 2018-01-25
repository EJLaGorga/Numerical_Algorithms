%A = zeros(19,2);
%for c = 1:24
%  A(c,1) = c/20;
%  A(c,2) = rosenbrock_linesearch_armijo(c/20)
%end

%plot(A(:,1), A(:,2))
%title('Convergence wrt \rho in Linesearch w/ Armijo')
%xlabel('\rho')
%ylabel('Iterations before ||-\Deltaf(x)|| <= 10^{-5}')

%linesearch_armijo([1.2,1.2]', @eval_rosenbrock, @grad_rosenbrock)
%linesearch_armijo([-1.2,1]', @eval_rosenbrock, @grad_rosenbrock)
%linesearch_armijo([1.2,1.2]', @eval_func, @grad_func)
linesearch_armijo([0,0]', @eval_func, @grad_func)

function k = linesearch_armijo(x, f_eval, f_grad)
  % parameters for linesearch with armijo rule
  fx = f_eval(x);           %evaluate f(x) at initial x
  P = f_grad(x);            %evaluate f'(x) at initial x
  p = .1;                   %random scalar in (0,1)
  c = .1;                   %random scalar in (0,1)
  k = 0;                    %counter

  % while the gradient is still sufficiently meaningful
  % which is equivelant to x sufficiently close to minimizer
  while  norm(P) > 1e-6      
    a = 1;
    fx_k = f_eval(x - a*P);
    while fx_k > fx + c*a*(P' * -P)
        a = p*a;
        fx_k = f_eval(x - a*P);
    end
    x = x - a*P;
    fx = fx_k;
    P = f_grad(x);
    k = k + 1;
  end
end

% evaluates the Rosenbrock function at a point in R^2
% takes 2d vector, returns scalar value
function y = eval_rosenbrock(x)
  y = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
end

% evaluates the gradient of the Rosenbrock function 
% at a point in R^2
% takes 2d vector, returns 2d vector
function y = grad_rosenbrock(x)
  y(1) = 400*x(1)^3 - 400*x(1)*x(2) + 2*x(1) - 2;
  y(2) = -200*(x(1)^2 - x(2));
  y = y';
end

% evaluates f(x) = sqrt(xq^2 + 1) + sqrt(xq^2 + 1)
% at a point in R^2, takes 2d vector, returns scalar value
function y = eval_func(x)
  %y = sqrt(x(1)^2 + 1) + sqrt(x(2)^2 + 1);
  y = 0;
  for i = 1:length(b)
    y = y + (1/2)*x'*A{i}*x - b{i}'*x;
  end
  y = y/length(b);
end

% evaluates the gradient of 
% f(x) = sqrt(xq^2 + 1) + sqrt(xq^2 + 1) at a point in R^2
% takes 2d vector, returns 2d vector
function y = grad_func(x)
  %y(1) = x(1)/sqrt(x(1)^2 + 1);
  %y(2) = x(2)/sqrt(x(2)^2 + 1);
  %y = y';y = zeros(length(x),1);
  for i = 1:length(b)
    y = y + (1/2)*(A{i}+A{i}')*x - b{i};
  end
  y = y/length(b);
end



