% The boundary points.
a = -1;
b = 1;

% The parameters of the initial condition.
c = 2;
d = 0.1;
T = 0.4;

% The initial condition.
ga = @(t)(c - tanh((a + 0.5 - c .* t) / (2 * e)));
gb = @(t)(c - tanh((b + 0.5 - c .* t) / (2 * e)));
u0 = @(x)(c-tanh((x + 0.5) / (2 * d)));

% The time initialization and its stepsize.
t = 0;
dts = 0.001;

% The number of nodes in the interval.
P = [41, 81, 161, 321, 641];
L = length(P);

% The error vector and the step size vector needed to compute the convergence rate.
e = zeros(1, L);
dxs = zeros(1, L);

%{
Step 1: Solve the equation for dU/dt by multiplying by M^{-1} from the left.
Step 2: The initial value of dU/dt is known as we know the initial condition of U.
Step 3: Use this step together with the standard 4th order explicit Runge-Kutta method
        in order to find the time deriative in each time step. THIS GIVES US WHAT?
%}

for i = 1:L
    % Setting the current problem size.
    N = P(i);
    
    % Discretizing the spatial domain.
    x = linspace(a, b, N);
    dx = x(2);
    dxs(i) = dx;
    
    dt = dts * dx;
    t = 0:dt:0.4;
    
    % Creating the solution matrix.
    u = zeros(length(x), length(t));
    
    % Defining the intial values.
    v0 = zeros(length(x), 1);
    for j = 1:length(x)
        v0(j) = u0(x(j));
    end
    
    % Adding the initial values.
    u(:, 1) = v0;
    
    % Construct the mass matrix.
    M = massMatrixAssembler(x);
    
    % Construct the stiffness matrix.
    A = advectionMatrixAssembler(x);

    % Construct the dispersion matrix.
    S = diffusionMatrixAssembler(x);
    
%     while t < T
%         t = t + dt;
%     end
end

function M = massMatrixAssembler(x)
K = length(x) - 1;
M = zeros(K + 1, K + 1);

for i = 1:K
    h = x(i+1) - x(i);
    M(i, i) = M(i, i) + h / 3;
    M(i, i+1) = M(i, i+1) + h / 6;
    M(i+1, i) = M(i+1, i) + h / 6;
    M(i+1, i+1) = M(i+1, i+1) + h / 3;
end

end

function A = advectionMatrixAssembler(x)
K = length(x) - 1;
S = zeros(K+1, K+1);

for i = 1:K
    h = x(i+1) - x(i);
    n = [i i+1];
    A(n, n) = A(n, n)
end

end

function S = diffusionMatrixAssembler(x)
K = length(x) - 1;
S = zeros(K + 1, K + 1);

for i = 1:K
    h = x(i+1) - x(i);
    n = [i i+1];
    S(n, n) = S(n, n) + (1 / h) * [1 -1; -1 1]; % Be careful about dividing by h.
end

S(1, :) = 0;
S(1, 1) = 1;
S(end, :) = 0;
S(end, end) = 1;
    
end