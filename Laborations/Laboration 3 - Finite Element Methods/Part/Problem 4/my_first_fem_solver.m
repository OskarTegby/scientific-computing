function my_first_fem_solver()
a = 0; % The left end point of interval
b = 1; % Right
N = [2, 4, 16, 256];
hv = [1/2, 1/4, 1/16, 1/256]; % The mesh size

for i=1:length(N)
    h = hv(i);
    x = a:h:b; % The node coords
    A1 = my_stiffness_matrix_assembler1(x);
    A2 = my_stiffness_matrix_assembler2(x);
    B1 = my_load_vector_assembler1(x);
    B2 = my_load_vector_assembler2(x);
    % This is the ratio between the largest and smallest values.
    fprintf("The condition number of A1 is %g.\n", cond(A1));
    fprintf("The condition number of A2 is %g.\n", cond(A2));
    xi1 = A1 \ B1; % Solving the system of equations
    xi2 = A2 \ B2;
    plot(x, xi1); % Plotting the solution
    hold on;
    plot(x, xi2); % What's the difference?
    hold off;
    pause(3);
end
close all;
end