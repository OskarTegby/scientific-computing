function A = my_stiffness_matrix_assembler(x)
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

A(1, 1) = 1.e+6;
A(N+1, N+1) = 1.e+6;