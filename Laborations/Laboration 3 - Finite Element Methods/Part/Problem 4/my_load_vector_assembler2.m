function B = my_load_vector_assembler2(x)
%
% Returns the assembled load vector b.
% Input is a vector x of node coords.
%
N = length(x) - 1;
B = zeros(N+1, 1);

for i = 1:N
    h = x(i+1) - x(i);
    n = [i i+1];
    B(n) = B(n) + [f(x(i)); f(x(i+1))] * h/2;
end

B(1) = 4;
B(end) = 7;