%A = zeros(19,2);
%for c = 1:24
%  A(c,1) = c/20;
%  A(c,2) = rosenbrock_linesearch_armijo(c/20)
%end

%plot(A(:,1), A(:,2))
%title('Convergence wrt \rho in Linesearch w/ Armijo')
%xlabel('\rho')
%ylabel('Iterations before ||-\Deltaf(x)|| <= 10^{-5}')

gradient_descent([1.2,1.2]', @eval_func, @grad_func)
gradient_descent([-1.2,1]', @eval_func, @grad_func)

function k = gradient_descent(x, f_eval, f_grad)
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
  
  function z = eval_residual(y, x)
    for i = 1:size(A)
        z = z + w(i)*((norm(y-a(i))^2)/(norm(x-a(i))));
    end
  end

  

end
