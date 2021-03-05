% Matrix size.
N = 10;

% Construct a matrix.
A = diag(repmat(2,1,N), 0) + diag(repmat(-1,1,N-1), -1) + diag(repmat(-1,1,N-1), 1);

% Right-hand-side vector.
b = 5*rand(N,1);

%
% Jacobi matrices
%

L = tril(A,-1); % Lower triangular (Jacobi, not LU)
U = triu(A,1); % Upper triangular (J, not LU)

D = diag(A).*eye(size(A));
Dinv = 1./diag(A).*eye(size(A)); % Inverse of a diagonal matrix. Very cheap.

%
% The iterative framework.
%

x0 = zeros(N,1); % Initial guess.
tolerance = 1e-4;

i=0; % An iterator.
while (norm(A*x0 - b,inf) > tolerance)
    x_next = Dinv * (-(L+U)*x0 + b); % Compute x^(i+1) = ....
    x0 = x_next; % Update x0.
    
    disp(['iteration: ' num2str(i) ' :::: max residual: ' num2str(norm(A*x0 - b,inf))])
    i=i+1    
end

% % x_0 is now the solution.
% A*x0 - b

    