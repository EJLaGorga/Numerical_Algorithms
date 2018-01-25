function [Err] = rand_kaczmarz(A,b)
% Implementation of random Kaczmarz algorithm
% Input paramerters: A overdetermined matrix
%                    b nonhomogenous solution

%p = sqrt(sum(A.^2,2))/norm(sqrt(sum(A.^2,2)),1);
p = sum(A.^2, 2)/norm(sum(A.^2, 2),1);
pd = makedist('Multinomial', 'probabilities', p);

n = size(A,2);
x = zeros(n,1);
xs = A\b;
Err = [];
while norm(xs - x) > 1e-1
    i = random(pd);
    x = x+(((b(i,1)-A(i,:)*x)/(norm(A(i,:))^2))*A(i,:)');
    Err = [Err; norm(xs - x)];
end


