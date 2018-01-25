trust_region_w_dogleg([1.2,1.2]', @eval_rosenbrock, @grad_rosenbrock, @hess_rosenbrock)
trust_region_w_dogleg([-1.2,1]', @eval_rosenbrock, @grad_rosenbrock, @hess_rosenbrock)


function k = trust_region_w_dogleg(xk, eval_f, grad_f, hess_f)
    TRmax = 2;  %maximum initial trust region radius
    TRk = 1;    %trust radius
    gmax = 1/2; 
    gamma = 0;
    k = 0;
    
    while norm(xk - [1,1]) > 1e-6
        E = eval_f(xk);
        G = grad_f(xk);
        H = hess_f(xk);
        
        t = fminbnd(@dogleg, 0, 2);
        [~, pk] = dogleg(t);
        rho = (E - eval_f(xk + pk))/(-1*G'*pk - (1/2)*pk'*H*pk);

        if rho < gmax
            TRk = gmax*TRk;
        elseif rho > (1-gmax) && norm(pk) == TRk
            TRk = min(2*TRk, TRmax);
        end   
        
        if rho > gamma
            xk = xk + pk;
        end
        k = k+1;        
    end

    function [y, p] = dogleg(t)
        PkSD = -1 * ((G' * G)/(G' * H * G)) * G;
        % thigh...
        if t <= 1
            p = t*PkSD;
            if norm(p) > TRk
                p = (TRk/norm(p))*p;
            end
        % calf :)
        elseif 1 < t
            PkQN = -1 * inv(H) * G;
            p = PkSD + (t-1)*(PkQN - PkSD);
            if norm(p) > TRk
                p = (TRk/norm(p))*p;
            end
        end
        % evaluation of m()
        y = E + G'*p + (1/2)*p'*H*p;
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

% evaluates the hessian of the Rosenbrock function 
% at a point in R^2
% takes 2d vector, returns 2x2 array
function H = hess_rosenbrock(x)
  H(1,1) = -400*(x(2) - x(1)^2) + 800*x(1)^2 + 2;
  H(1,2) = -400*x(1);
  H(2,1) = -400*x(1);
  H(2,2) = 200;
end