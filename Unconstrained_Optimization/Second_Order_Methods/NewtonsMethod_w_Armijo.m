newtons_method_w_armijo([1.2,1.2]', @eval_func, @grad_func, @hess_func)
newtons_method_w_armijo([-1.2,1]', @eval_func, @grad_func, @hess_func)
newtons_method_w_armijo([1.2,1.2]', @eval_rosenbrock, @grad_rosenbrock, @hess_rosenbrock)
newtons_method_w_armijo([-1.2,1]', @eval_rosenbrock, @grad_rosenbrock, @hess_rosenbrock)

function k = newtons_method_w_armijo(x, eval_f, grad_f, hess_f)
    p = x - [0,0]';
    c = 0.1;
    a = 0.1;
    k = 0;
    while norm(p) > 1e-6
        R  = chol(hess_f(x));
        e  = R'\(-grad_f(x));
        p = R\e;
        while eval_f(x + a*p) > eval_f(x) + c*a*grad_f(x)'*p
            a = p*a;
        end
        x =  x + a*p;
        k = k + 1;
    end
end


% evaluates f(x) = sqrt(xq^2 + 1) + sqrt(xq^2 + 1)
% at a point in R^2, takes 2d vector, returns scalar
function y = eval_func(x)
  y = sqrt(x(1)^2 + 1) + sqrt(x(2)^2 + 1);
end

% evaluates the gradient of 
%f(x) = sqrt(xq^2 + 1) + sqrt(xq^2 + 1) at a point in R^2
% takes 2d vector, returns 2d vector
function y = grad_func(x)
  y(1) = x(1)/sqrt(x(1)^2 + 1);
  y(2) = x(2)/sqrt(x(2)^2 + 1);
  y = y';
end

% evaluates the hessian of 
% f(x) = sqrt(xq^2 + 1) + sqrt(xq^2 + 1) at a point in R^2
% takes 2d vector, returns 2x2 matrix
function H = hess_func(x)
  H(1,1) = 1/sqrt(x(1)^2 + 1) - (x(1)^2)/((x(1)^2 + 1)^(3/2));
  H(1,2) = 0;
  H(2,1) = 0;
  H(2,2) = 1/sqrt(x(2)^2 + 1) - (x(2)^2)/((x(2)^2 + 1)^(3/2));
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

% evaluates the hessian of the Rosenbrock function 
% at a point in R^2
% takes 2d vector, returns 2x2 array
function H = hess_rosenbrock(x)
  H(1,1) = -400*(x(2) - x(1)^2) + 800*x(1)^2 + 2;
  H(1,2) = -400*x(1);
  H(2,1) = -400*x(1);
  H(2,2) = 200;
end

