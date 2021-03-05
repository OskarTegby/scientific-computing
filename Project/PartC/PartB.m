%% The solution of Part B
while true
part = 0;
toggle = 1;
while part ~= 1 && part ~= 2 && part ~= 4
    part = input("Please enter which part to run (1, 2, or 4) or 0 to quit: ");
    if part == 0
        toggle = 0;
        break;
    end
end
if toggle == 0
    break;
end

check = input("Plot the solutions (0/1). (Default 1): ");
if isempty(check)
    check = 1;
end

% The different problem sizes.
M = [40, 80, 160, 320, 640];
L = length(M);

% The space starts at a and ends at b.
a = 0;
b = 1;

% The time starts at 0 and ends at T. 
T = 1;
dts = [0.01 0.001];
dxs = zeros(1, L);

% The error vector.
e = zeros(1, L);

for i = 1:L
    % Setting the current problem size.
    N = M(i);
    
    % Discretizing the spatial domain.
    x = linspace(a, b, N);
    dx = x(2);
    dxs(i) = dx;
    
	% Discretizing the temporal domain.
    dt = dts(1) * dx;
    t = 0:dt:1;

    % Defining the piecewise function.
    syms X;
    if part == 1
        u0 = piecewise(abs(2 * X - 0.3) <= 0.25, exp(-300 * (2 * X - 0.3).^2), 0);
    elseif part == 2 || part == 4
        u0 = piecewise(abs(2 * X - 0.3) <= 0.25, exp(-300 * (2 * X - 0.3).^2),...
                       abs(2 * X - 0.9) <= 0.2, 1,...
                       abs(2 * X - 1.6) <= 0.2, sqrt(1-((2 * X - 1.6)/0.2)^2), 0);
    end
    
    % Defining the space-discretization matrix.
    if part == 1 || part == 2
        A = 1 / (2 * dx) * spdiags(repmat([-1, 1, -1, 1], N, 1), [-(N-1) -1, 1, (N-1)], N, N);
    elseif part == 4
        A = (-1 / dx) * spdiags(repmat([1, -1, 0, -1], N, 1), [0, -1, 1, N-1], N, N);
    end
        
    % Creating the solution matrix.
    u = zeros(length(x), length(t));
    
    % Adding the initial values.
    u(:, 1) = subs(u0, X, x(1:end));
    
    % Performing Runge-Kutta 4 on the solution.
    for j = 2:length(t)
        k1 = A * u(:, j-1);
        k2 = A * u(:, j-1) + dt * k1 / 2;
        k3 = A * u(:, j-1) + dt * k2 / 2;
        k4 = A * u(:, j-1) + dt * k3;
        u(:, j) = u(:, j-1) + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    end
        
    % Computing the error.
    e(i) = norm((u(:, 1) - u(:,length(t)))/sqrt(N));
%     fprintf("The error for the problem size %d is %.2g.\n", N, e(i));
    
    if part == 1 || part == 2
        s = 100;
    elseif part == 4
        s = 200;
    end
    figure;
    if check
        for k = 1:s:length(t)
            hold off;
            plot(x, u(:,k), '-k', 'LineWidth', 2);
            title("Simulation for N = "+num2str(N)+".");
            xlabel("x");
            ylabel("u(x)");
            axis([x(1), x(end), -1.5, 1.5]);
            drawnow;
        end
        close all;
    end
end

close all;

M = flip(M);
figure;
loglog(1./M, e, 'LineWidth', 2);
hold on
xlabel("The logarithm of the problem size, log(N).");
ylabel("The logarithm of the error in the L2 norm, |u-u_h|");
title("The loglog plot of the error");
hold off
% exportgraphics(gca, "Part B4 - loglog plot.png");
% exportgraphics(gca, "Part B4 - loglog plot.pdf");
fprintf("Press enter to continue.\n");
pause;
close all;

rate = zeros(1, L-1);
for i = 2:L
    % Derivation of the convergence rate:
    % e1 / e2 = (h1 / h2)^rate
    % log(e1 / e2) = rate*log(h1 / h2)
    % rate = log(e1 / e2) / log(h1 / h2)
    rate(i-1) = (log(e(i)) - log(e(i - 1))) / (log(dxs(i)) - log(dxs(i - 1)));
end

plot(2:L, rate);
xlabel("Index");
ylabel("Convergence rate");
title("The convergence rate for each step size");
fprintf("Press enter to continue.\n");
pause;
close all;
% exportgraphics(gca, "Part B4 - Convergence rate.png");
% exportgraphics(gca, "Part B4 - Convergence rate.pdf");
end