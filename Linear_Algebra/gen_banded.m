function M = gen_banded(n,b)
%generate a matrix of size n with semibandwidth of b
M = sign(conv2(eye(n),ones(b+1),'same'));
for i = 1 : n
	for j = 1 : n
		if M(i,j) ~= 0
			M(i,j) = rand();
		end;
	end;
end;
