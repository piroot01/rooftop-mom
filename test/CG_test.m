clc; clear; format long;

%--- Problem Setup ---
n = 100;
A = rand(n);
B_vec = rand(n, 1);
x = zeros(n, 1);

maxIter = 100000;
tol = 1e-5;

%--- Precompute Transpose ---
A_trans = A.';

%--- Initial Residuals ---
r = B_vec - A * x;
res_sq = r' * r;  % Squared residual norm
tol_sq = tol^2;

r_tld = r;  % Shadow residual
p = r;
p_tld = r_tld;

rho = r_tld' * r;

%fprintf('Iteration %d: residual norm = %e\n', 0, sqrt(res_sq));

%--- Main BiCG Loop ---
iter = 0;
while (res_sq > tol_sq) && (iter < maxIter)
    Ap = A * p;
    Atp_tld = A_trans * p_tld;  % Use precomputed transpose
    
    denom = p_tld' * Ap;
    if abs(denom) < eps
        error('Breakdown: denominator too small.');
    end
    
    alpha = rho / denom;
    
    % Update solution and residuals
    x = x + alpha * p;
    r = r - alpha * Ap;
    r_tld = r_tld - alpha * Atp_tld;
    
    % Compute new residual and rho
    res_sq_new = r' * r;
    rho_new = r_tld' * r;
    sqrt_res = sqrt(res_sq_new);
    
    % Check convergence
    if res_sq_new < tol_sq
        iter = iter + 1;
        fprintf('Iteration %d: residual norm = %e\n', iter, sqrt_res);
        break;
    end
    
    % Update search directions
    beta = rho_new / rho;
    p = r + beta * p;
    p_tld = r_tld + beta * p_tld;
    
    % Prepare for next iteration
    rho = rho_new;
    res_sq = res_sq_new;
    iter = iter + 1;
    %fprintf('Iteration %d: residual norm = %e\n', iter, sqrt_res);
end


function y = applyOperator(x, flag, A)
    if strcmp(flag,'notransp') % Compute A*x
        y = A * x;
    elseif strcmp(flag,'transp') % Compute A'*x
        y = A.' * x;
    end
end

tic;

x = bicg(@(x, flag) applyOperator(x, flag, A), B_vec, tol, maxIter);

toc;

%--- Validation ---
x_exact = A \ B_vec;
fprintf('Exact solution norm: %e\n', norm(x_exact));
fprintf('BiCG solution norm: %e\n', norm(x));
