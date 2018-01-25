function y = evaluate_f(A, b, x)
y = 0;
for i = 1:length(b)
    y = y + (1/2)*x'*A{i}*x - b{i}'*x;
end
y = y/length(b);
end