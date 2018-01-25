function y = gradient_f(A, b, x)
y = zeros(length(x),1);
for i = 1:length(b)
    y = y + (1/2)*(A{i}+A{i}')*x - b{i};
end
y = y/length(b);
end