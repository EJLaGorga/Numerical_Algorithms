function [L,A]=lu_nopivot(A)
% performs gausian elimination on matrix resulting 
% in an upper triangular matrix A, and a lower 
% triangular matix L

n = size(A,1);
L = eye(n);
for i = 1 : n
	for j = i+1 : n
		L(j,i) = A(j,i)/A(i,i);
		for k = i : n
			A(j,k) = A(j,k) - A(i,k) * L(j,i);
		end;
	end;
end;	
