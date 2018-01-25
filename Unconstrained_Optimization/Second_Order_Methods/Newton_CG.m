for i = 2:6
    xk = -1*ones(2^i, 1);
    [itters(i-1, 1), ~] = Newton_ConjugateGradient(xk);
    itters(i-1, 2) = 2^i;
end

figure
plot(itters(:,2), itters(:,1))
title('Newton\_CG Convergence Results \gamma_{k} = min(.5, ||\DeltaF(x)||)')
%title('Newton\_CG Convergence Results \gamma_{k} = .1')
xlabel('dimensions of objective functions domain')
ylabel('iterations to convergence')

xk = -1*ones(2^3, 1);
[~, dxk] = Newton_ConjugateGradient(xk);
dxk'

function [k, XK] = Newton_ConjugateGradient(xk)
% f:Rn->R [evaluates function at point]
evaluate = @(x) sum(100*(x(2:2:end)-x(1:2:end).^2).^2+(1-x(1:2:end)).^2);
% f:Rn->Rn [evaluates gradient at point]
dxodd =  @(x) 400*(x(1:2:end).^3 - x(2:2:end).*x(1:2:end))+2*x(1:2:end)-2;
dxeven = @(x) 200*(x(2:2:end)-x(1:2:end).^2);
gradient = @(x) reshape([dxodd(x)'; dxeven(x)'], 1, [])';
% f:Rn->Rnxn [evaluates hessian at point]
dxoo = @(x) 1200*x(1:2:end).^2-400*x(2:2:end)+2;
dxee = @(x) 200*ones(length(x)/2, 1);
diagonal = @(x) reshape([dxoo(x)'; dxee(x)'], 1, [])';
dxe = @(x) -400*x(1:2:end);
dxo = @(x) zeros(length(x)/2, 1);
remf = @(x) x(2:end);
tridiag = @(x) remf(reshape([dxo(x)'; dxe(x)'], 1, []))';
hessian = @(x) diag(diagonal(x),0) + diag(tridiag(x), -1) + diag(tridiag(x), 1);

% Newtons algorithm (solve system, update xk using armijo)
    XK = [];
    p = .1;
    c = .5;
    k = 0;
    while norm(gradient(xk)) > 1e-6
        Gk = gradient(xk);
        gma = min(0.5, sqrt(norm(Gk)));
        %gma = .1;
        Hk = hessian(xk);
        Pk = CG_Inexact_Newton(Hk, Gk);
        
        fk = evaluate(xk);
        a = 1;
        while evaluate(xk+a*Pk) > fk + c*a*(Gk'*Pk)
            a = p*a;
        end
        
        XK = [XK, norm(xk + a*Pk - ones(length(xk)))/norm(xk - ones(length(xk)))];
        
        xk = xk + a*Pk;
        k = k + 1;
    end

    function Pk = CG_Inexact_Newton(Hk, Gk)
        % Initialize residual vector, and conjugate direction
        zj = zeros(length(Gk), 1);
        rj = Gk;
        dj = -rj;
        for j = 1:length(Gk)
            if dj'*Hk*dj < 0
                if j == 1
                    Pk = -Gk; return;
                else
                    Pk = zj; return;
                end
            end
            aj = (rj'*rj)/(dj'*Hk*dj);
            zj = zj + aj*dj;
            rj1 = rj + aj*Hk*dj;
            if norm(rj1) <= gma*norm(Gk)
                Pk = zj; return;
            end
            dj = -rj1 + ((rj1'*rj1)/(rj'*rj))*dj;
            rj = rj1;
        end
    end
end