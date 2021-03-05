function A = my_stiffness_matrix_assembler(x)
%
% Returns the assembled stiffness matrix A.
% Input is a vector x of node coords.
%
N = length(x) - 1; % The number of elements
A = zeros(N+1, N+1); % Initializing the stiffness matrix to zero
for i = 1:N % Looping over the elements
    h = x(i+1) - x(i); % The element lengths
    n = [i i+1]; % The nodes
    A(n, n) = A(n, n) + [1 -1; -1 1]; % Assembling the element stiffness
end

A(1, :) = 0;
A(1, 1) = 1;
A(N+1, :) = 0;
A(N+1, N+1) = 1;