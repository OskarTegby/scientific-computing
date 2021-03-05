function my_first_fem_solver()
a = 0;
b = 1;
h = 1/2;

x = a:h:b;
A = my_stiffness_matrix_assembler(x);
B = my_load_vector_assembler(x);
xi = A \ B;
plot(x, xi);
pause(3);

end