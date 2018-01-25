A = zeros(4,1);
A(1,1) = conjugate_gradient(zeros(5,1), ones(5,1), hilb(5));
A(2,1) = conjugate_gradient(zeros(8,1), ones(8,1), hilb(8));
A(3,1) = conjugate_gradient(zeros(12,1), ones(12,1), hilb(12));
A(4,1) = conjugate_gradient(zeros(20,1), ones(20,1), hilb(20));
A

function k = conjugate_gradient(x, b, A)
    % Initialized count, residual vector, and conjugate direction
    k = 0;
    r = b - A*x;
    p = r;
    while norm(r) > 1e-6
        % Initialize ak and calculate xk at each itteration
        a = (r'*r)/(p'*A*p);
        x = x + a*p;
        % Calculate rk+1 and initialize Bk at each itteration
        r_1 = r - a*A*p;
        B = (r_1'*r_1)/(r'*r);
        r = r_1;
        % Calculate pk with Bk 
        p = r + B*p;
        % Itterate count
        k = k + 1;
    end
end
