while true
%% Specifying which parts to run.
run = input("Please enter which part to run (1, 2, or 3) or 0 to quit: ");
if ismember(run, [1, 2, 3]) == false
    return;
end
    
check = input("Check against divergent executions. (0/1). (Default 1): ");
if isempty(check)
    check = 1;
end
verbose = input("Label the outputs and exceptions (0/1). (Default 0): ");
if isempty(verbose)
    verbose = 0;
end

%% Part A1
if run == 1
    A = [10 -1 2 0; -1 11 -1 3; 2 -1 10 -1; 0 3 -1 8];
    b = [6; 25; -11; 15];
    TOL = 0.001;
    
    fprintf("\n============================================\n");
    fprintf("\n             Performing A1\n");
    fprintf("\n============================================\n");
    calls(A, b, TOL, run, check, verbose);
    fprintf("\n");
end

%% Part A2
if run == 2
    sizes = [100, 500, 1000];
    weights = [1, 5, 10, 100];
    TOL = 0.00001;
    safe = 0;

    fprintf("\n============================================\n");
    fprintf("               Performing A2");
    fprintf("\n============================================\n");
    for i = 1:length(sizes)
        for j = 1:length(weights)
            N = sizes(i);
            w = weights(j);
            fprintf("Running the solvers with N = %g with w = %g.\n", N, w);
            A = rand(N) + diag(w * ones(N, 1));
            
            % Checks that the generated matrix qualifies is safe == 1.
            if safe && issymmetric(A) && all(eig(A) > 0)
                run = -1;
            elseif safe == 0
                run = -1;
            end
            b = rand(N, 1);
            calls(A, b, TOL, run, check, verbose);
            fprintf("============================================\n");
        end
    end
end

%% Part A3
if run == 3
    alpha = [1, 0.1, 0.001, 0.00001];
    TOL = 0.00001;
    N = 100000;
    b = rand(N, 1);

    fprintf("\n============================================\n");
    fprintf("               Performing A3");
    fprintf("\n============================================\n");
    for i = 1:length(alpha)
        fprintf("Running the solvers with alpha = %g.\n", alpha(i));
        A = spdiags(repmat([-1,2+alpha(i),-1], N ,1),[1,0,-1], N, N);
        calls(A, b, TOL, run, check, verbose);
        fprintf("============================================\n");
    end
end

end


%% This is the function which calls all the solvers.
function calls(A, b, TOL, run, check, verbose)

A = full(A);
if verbose == 0
    results = zeros(6, 2);
end

tic;
A = sparse(A);
[x2, i] = jacobi(A, b, TOL, check); %#ok<ASGLU>

if verbose
    if isnan(i)
        fprintf("The solution diverged using Jacobi iteration.\n");
    else
        fprintf("Jacobi iteration took %g seconds and %g iterations.\n", toc, i);
    end
else
    if isnan(i)
        results(1, 1) = NaN;
        results(1, 2) = NaN;
    else
        results(1 ,1) = i;
        results(1, 2) = toc;
    end
end

tic;
[x3, i] = gs(A, b, TOL, check); %#ok<ASGLU>

if verbose
    if isnan(i)
        fprintf("The solution diverged using the Gauss-Seidel method.\n");
    else
        fprintf("The Gauss-Seidel method took %g seconds and %g iterations.\n", toc, i);
    end
else
    if isnan(i)
        results(2, 1) = NaN;
        results(2, 2) = NaN;
    else
        results(2, 1) = i;
        results(2, 2) = toc;
    end
end

if run ~= 2 || check ~= 1
    tic;
    [x4, i] = cg(A, b, TOL, check); %#ok<ASGLU>
    
    if verbose
        if isnan(i)
            fprintf("The solution diverged using the Conjugate-Gradient method.\n");
        else
            fprintf("The Conjugate-Gradient method took %g seconds and %g iterations.\n", toc, i);
        end
    else
        if isnan(i)
            results(3, 1) = NaN;
            results(3, 2) = NaN;
        else
            results(3, 1) = i;
            results(3, 2) = toc;
        end
    end
else
    if verbose
        fprintf("The Conjugate-Gradient method was skipped as the matrix did not qualify.\n");
    else
        results(3, 1) = NaN;
        results(3, 2) = NaN;
    end
end

tic;
x5 = myLU(A, b); %#ok<NASGU>

if verbose
    fprintf("My LU solver took %g seconds.\n", toc);
else
    results(4, 1) = 1;
    results(4, 2) = toc;
end

tic;
[L, U] = lu(A);
y = L \ b;
x6 = U \ y; %#ok<NASGU>

if verbose
    fprintf("MATLAB's LU solver took %g seconds.\n", toc);
else
    results(5, 1) = 1;
    results(5, 2) = toc;
end

tic; 
x1 = A \ b; %#ok<NASGU>

if verbose
    fprintf("The backslash operator took %g seconds.\n", toc);
else
    results(6, 1) = 1;
    results(6, 2) = toc;
end

if verbose == 0
    results %#ok<NOPRT>
end

end

%% This performs Jacobi iteration.
function [x, i] = jacobi(A, b, TOL, check)
    s = size(A,1);
    d = (1./diag(A))';
    Dinv = sparse(1:s, 1:s, d, s, s);
    L = tril(A, -1);
    U = triu(A, 1);
    LU = L + U;
    M = -Dinv * LU;
    c = Dinv * b;
    x = zeros(s, 1); % The initial guess
    i = 0; % The number of iterations

    disp(norm(full(-Dinv) * full(LU)));
    
    if (ismember(1, all(norm(-full(Dinv) * full(LU)) < 1)))
        fprintf("The method converges.\n");
    end
    
    Err = 1;
    while Err > TOL
        Err = norm(A*x - b);
        x =  M * x + c;
        i = i + 1;
        
        if check % The check parameter disables the divergence checks for performance gain.
            if ismember(1, isnan(x)) || ismember(-Inf, x) || ismember(Inf, x)
                i = NaN;
                break;
            end
        end
    end
end

%% This performs the Gauss-Seidel method.
function [x, i] = gs(A, b, TOL, check)
    L = tril(A, -1);
    U = triu(A, 1);
    D = diag(diag(A));
    DL = D+L;
    DLinv = inv(DL);
    DLb = DL \ b;
    x = zeros(length(A), 1);
    i = 0;
    
    disp(norm(-full(DLinv) * full(U)));
    
    if (ismember(1, all(norm(-full(DLinv) * full(U)) < 1)))
        fprintf("The method converges.\n");
    end
    
    while norm(A*x - b) > TOL
        x = DLb + DL \ (-U * x);
        i = i + 1;
        
        if check
            if ismember(1, isnan(x)) || ismember(-Inf, x) || ismember(Inf, x)
                i = NaN;
                break;
            end
        end
    end 
end

%% This performs the conjugate gradient method.
function [x, i] = cg(A, b, TOL, check)
	x = zeros(length(A), 1);
    r = b - A * x;
    rho = r' * r;
    k = 0;
    i = 0;
    timeout = 10000;
    
    while (sqrt(rho) > TOL)
        k = k + 1;
        if k == 1
            P = r;
        else
            beta = rho / rho2;
            P = r + beta * P;
        end
        w = A * P;
        alpha = (P' * r) / (P' * w);
        x = x + alpha * P;
        r = r - alpha * w;
        rho2 = rho; 
        rho = r' * r;
        i = i + 1;
        
        if check || i >= timeout
            if (ismember(1, isnan(x)) || ismember(-Inf, x) || ismember(Inf, x)) || i >= timeout
                i = NaN;
                break;
            end
        end
    end
end

%% This is my own implementation of the LU decomposition algorithm.
function x = myLU(A, b)
    A = full(A);
    U = A;
    N = length(A);
    L = diag(ones(N,1));
    
    for k = 1:N-1
        for i = k+1:N
            if U(i,k) ~= 0
                factor = U(i, k) / U(k, k);
                L(i, k) = factor;
                for j = 1:N
                    U(i, j) = U(i, j) - factor * U(k, j);
                end
            end
        end
    end
    d = L \ b;
    x = U \ d;
end