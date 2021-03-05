clearvars;
close all;
%
% Use FDM to solve u_t = u_xx, x in (0,1), u(0,t)=0, u(1,t)=0, u(x,0) = sin(2*pi*x)
%

% Initial condition.
u_init = @(x) sin(2*pi*x);

%
% Input parameters.
%
x0 = 0; % End of the spatial domain.
xEnd = 1; % Beginning of the spatial domain.

tEnd = 20; % The final time.
dt = 0.001; % Time step.

u_bc = [0 0]; % Vector of the boundary conditions.

N = 30; % Choose the amount of points you want to spread over Omega.

%
% Discretize the space.
%

% (1) Discretize the spatial domain Omega.
X = linspace(x0,xEnd,N)'; % Equidistant points on [x0, xEnd].

% (2) Form the sparse matrix A.
h = X(2)-X(1); % Recover internodal distance from the set of points.
A = 1/h^2 * spdiags(repmat([-2, 1, 1], N-2, 1), [0, -1, 1], N-2, N-2);

%
% Now solve the system of ODEs u_t = A*u in time.
%

% (3) Discretize the time domain.
tt = 0:dt:tEnd;

% (4) Sample the initial condition.
u0 = u_init(X(2:end-1));

% (5) Apply Euler forward to solve.
u = zeros(length(X)-2, length(tt));
u(:,1) = u0;

for k = 2:length(tt)
    u(:,k) = u(:,k-1) + dt*A*u(:,k-1);
end


% (6) Draw the solution in time.
u = [u_bc(1)*ones(size(tt)); u; u_bc(2)*ones(size(tt))]; % Just append the BCs for every time step.

figure;
for k=1:length(tt)
    hold off;
    plot(X, u(:,k), '-k', 'LineWidth', 2);
    axis([X(1), X(end), -1, 1]);
    drawnow;
end
