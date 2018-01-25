function [L, A, P] = lu_partpivot(A)
n = size(A); 
L = eye(n); 
P = eye(n); 

for i = 1 : n
	[pivot m] = max(A(i:n, i));     
	m = m + i-1;
 	if m ~= i
		%swap rows
 		A([m,i],:) =  A([i,m], :);   
		P([m,i],:) =  P([i,m], :);   
		if i >= 2;          
			L([m,i], 1:i-1) =  L([i,m], 1:i-1);    
		end;
	end;

	for j = i+1 : n      
		%calculates gausian coefficient puts at index in l
		L(j, i) = A(j, i) / A(i, i);
		for k = i : n
			%calculates the the nth iteration of A
    			A(j, k) =  A(j, k) - L(j, i) * A(i, k);
  		end;
	end;
end;

