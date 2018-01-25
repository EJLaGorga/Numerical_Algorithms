function [Err] = kaczmarz(A, b)
% Implementation of Kaczmarz algorithm
% Input paramerters: A overdetermined matrix
%                    b nonhomogenous solution
m = size(A,1);
n = size(A,2);
x = zeros(n,1);
xs = A\b;
Err = [];
while norm(xs - x) > 1e-1
    i = mod(length(Err), m)+1;
    x = x+(((b(i,1)-A(i,:)*x)/(norm(A(i,:))^2))*A(i,:)');
    Err = [Err; norm(xs - x)];
end



