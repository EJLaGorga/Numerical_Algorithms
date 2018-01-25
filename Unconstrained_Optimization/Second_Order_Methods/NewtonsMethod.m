newtons_method([1.2,1.2]', @grad_rosenbrock, @hess_rosenbrock)
newtons_method([-1.2,1]', @grad_rosenbrock, @hess_rosenbrock)
newtons_method([1.2,1.2]', @grad_func, @hess_func)
newtons_method([-1.2,1]', @grad_func, @hess_func)

function k = newtons_method(x, grad_f, hess_f)
    dx = x - [1,1]';
    k = 0;
    while norm(dx) > 1e-6
        % "NEVER Invert a matrix" - Misha Kilmer, NLA 2017
        R  = chol(hess_f(x));
        c  = R'\(-grad_f(x));
        dx = R\c;
        x = dx + x;
        k = k + 1;
    end
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





