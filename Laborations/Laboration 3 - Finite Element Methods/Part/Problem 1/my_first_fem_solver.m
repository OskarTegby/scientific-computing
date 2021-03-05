function my_first_fem_solver()
a = 0; % The left end point of interval
b = 1; % right
N = [2, 4, 16, 256];
hv = [1/2, 1/4, 1/16, 1/256]; % The mesh size

for i=1:length(N)
    h = hv(i);
    x = a:h:b; % The node coords
    A = my_stiffness_matrix_assembler(x);
    B = my_load_vector_assembler(x);
    % This is the ratio between the largest and smallest values.
    fprintf("The condition number of A is %g.\n", cond(A));
%     fprintf("The condition number of B is %g.\n", cond(B));
    xi = A \ B; % Solving system of equations
    plot(x, xi); % Plotting solution
    hold on
    plot(x, x.*(1-x)); % Why is this so much smaller?
    hold off
    pause(3);
end
close all;
end