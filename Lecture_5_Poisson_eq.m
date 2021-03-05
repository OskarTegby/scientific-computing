clearvars;
close all;
%
% Use FDM to solve u_xx = -4*pi^2*sin(2*pi*x), x in (0,1), u(0)=0, u(1)=0
%

% Exact solution to this equation:
u_exact = @(x) sin(2*pi*x);

%
% Input parameters.
%
f = @(x) -4*pi^2*sin(2*pi*x); % The right-hand-side function (forcing function).
x0 = 0; % End of the domain.
xEnd = 1; % Beginning of the domain.
u_bc = [0 0]; % Vector of the boundary conditions.

N = 10; % Choose the amount of points you want to spread over Omega.

%
% Compute.
%

% (1) Discretize the domain Omega.
X = linspace(x0,xEnd,N)'; % Equidistant points on [x0, xEnd].

% (2) Sample the RHS forcing function.
b = f(X(2:end-1)); % 2:end-1, because we are not solving the system for known BCs.

% (3) Impose the known boundary conditions.
h = X(2)-X(1); % Recover internodal distance from the set of points.
b(2) = b(2) - u_bc(1)/h^2;
b(end) = b(end) - u_bc(2)/h^2;


% (4) Form the sparse matrix A.
A = 1/h^2 * spdiags(repmat([-2, 1, 1], N-2, 1), [0, -1, 1], N-2, N-2);

% (5) Solve.
u = A\b;

% (6) Plot.
u = [u_bc(1); u; u_bc(2)]; % Append the boundary conditions for visualization.
figure;
plot(X, u, 'xk'); hold on;
plot(X, u_exact(X), 'or');


% (7) Compare to true solution.
norm(u - u_exact(X), 2)
