function [Err] = biasedSGD(A,b)
% Implementation of SGD with partially biased sampling
% Input paramerters: A overdetermined matrix
%                    b nonhomogenous solution

n = length(b);
w = (1/2*n) + (1/2)*sum(A.^2,2)/sum(sum(A.^2));
pd = makedist('Multinomial', 'probabilities', w/norm(w,1));
c = .1;
Af = sum(sum(A.^2));

x = zeros(size(A,2),1);
xs = A\b;
Err = [];
while norm(xs - x) > 1e-1
    i = random(pd);
    x = x + 2*c*(b(i) - A(i,:)*x)/(Af/n + sum(A(i,:).^2))*A(i,:)';
    Err = [Err; norm(xs - x)];
end

end