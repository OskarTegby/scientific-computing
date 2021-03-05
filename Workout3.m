N = 4;
A = spdiags(repmat([1, -1, -1, 1], N, 1), [0, -1, 1, 0], N, N);

syms x;
b = zeros(N, 1);
for i=1:N
    b(i) = 40 * int((x - (i / 4)) * sin(x), 0, 1);
end
b
u = zeros(1, N);

for i=1:N
    u(i) = 10 * sin(i / 4) - 10 * sin(1) * i / 4;
end

plot(0.25:0.25:1, c);
hold on
plot(0.25:0.25:1, u);
hold off
xlabel("x");
ylabel("y");
legend("u_h", "u");
title("The finite element solution");