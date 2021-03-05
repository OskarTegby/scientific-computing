function A = my_stiffness_matrix_assembler2(x)
%
% Returns the assembled stiffness matrix A.
% Input is a vector x of node coords.
%
N = length(x) - 1;
A = zeros(N+1, N+1);
for i = 1:N
    h = x(i+1) - x(i);
    n = [i i+1];
    A(n, n) = A(n, n) + [1 -1; -1 1];
end

A(1, :) = 0;
A(1, 1) = 1;
A(N+1, :) = 0;
A(N+1, N+1) = 1;