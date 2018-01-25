
m = 2;
n = 4;
xk = -1*ones(n, 1);
lm_gd_armijo(xk, m, n)

% Newtons algorithm (solve system, update xk)
function k = lm_gd_armijo (xk, m, n)
    evaluate = @(x) sum(100*(x(2:2:end)-x(1:2:end).^2).^2+(1-x(1:2:end)).^2);

    dxodd =  @(x) 400*(x(1:2:end).^3 - x(2:2:end).*x(1:2:end))+2*x(1:2:end)-2;
    dxeven = @(x) 200*(x(2:2:end)-x(1:2:end).^2);
    gradient = @(x) reshape([dxodd(x)'; dxeven(x)'], 1, [])';
    
    memory = zeros(n,m,2);
    
    p = .1;
    c = .5;
    k = 0;
    while norm(gradient(xk)) > 1e-6
        s1 = memory(:,1,1);
        y1 = memory(:,1,2);
        if y1'*y1 == 0
            Hko = eye(length(xk));
        else
            Hko = ((s1'*y1)/(y1'*y1))*eye(length(xk));
        end
            
        Gk = gradient(xk);    
        Pk = lm_BFGS(Gk);
        
        fk = evaluate(xk);
        a = 1;
        while evaluate(xk+a*Pk) > fk + c*a*(Gk'*Pk)
            a = p*a;
        end
    
        sk = a*Pk;
        yk = gradient(xk + a*Pk) - gradient(xk);
    
        memory(:, 2:m,:) = memory(:,1:m-1,:);
        memory(:,1,1) = sk;
        memory(:,1,2) = yk;
    
        xk = xk + a*Pk;
        k = k + 1;
    end
    
    function Pk = lm_BFGS(qk)
        for i = 1:m
            si = memory(:,i,1);
            yi = memory(:,i,2);
            if yi'*si == 0; pi = 0; else; pi = 1/(yi'*si); end
            ai = pi*si'*qk;
            qk = qk - ai*yi;
        end
        rk = Hko*qk;
        for i = m:-1:1
            si = memory(:,i,1);
            yi = memory(:,i,2);
            if yi'*si == 0; pi = 0; else; pi = 1/(yi'*si); end
            ai = pi*si'*qk;
            Bi = pi*yi'*rk;
            rk = rk + si*(ai - Bi);
        end
        Pk = -rk;
    end
end






