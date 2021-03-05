function my_first_fem_solver()
a = 0;
b = 1;
N = [2, 4, 16, 256];
hv = [1/2, 1/4, 1/16, 1/256];
e = zeros(length(N));
rate = zeros(length(N));

for i = 1:length(N)
    h = hv(i);
    x = a:h:b;
    A = my_stiffness_matrix_assembler(x);
    B = my_load_vector_assembler(x);
    xi = A \ B; % Why does it increase so much?
    uh2 = xi' * A * xi;
    u2 = sum((1 - 2 * x).^2);
    e(i) = u2 - uh2;
    plot(x, xi);
    hold on
    plot(x, x.*(1-x));
    xlabel("x");
    ylabel("u(x)");
    legend("u_h", "u");
    title("Analytical and numerical solutions.");
    hold off
%     pause(5);
end

for i = 2:length(N)
    rate(i) = (log(hv(i - 1)) - log(hv(i))) / (log(e(i - 1)) - log(e(i)));
end

plot(1:length(N), rate);
xlabel("i");
ylabel("Convergence rate");
title("The convergence rate of the simulation.");
% pause(60);
plot(1:length(N), e);
xlabel("i");
ylabel("Error");
title("The error of the simulation.");

pause(60);
close all;
end