clc; clear; format long;

%--- Problem Setup: Construct an SPD matrix A ---
n = 100;                        % You can choose an appropriate size n
R = rand(n);                  % Generate a random matrix R
A = R' * R;                   % A is now symmetric positive definite
B_vec = rand(n,1);            % Right-hand side vector

%--- Define operator and initial guess ---
applyOperator = @(x) A * x;
J_vec = zeros(n,1);           % Initial guess (can be zeros or any other guess)

%--- CG Parameters ---
maxIter = 1000;
tol = 1e-3;

tic;

%--- Compute initial residual ---
r = B_vec - applyOperator(J_vec);
p = r;
res_norm = norm(r) / norm(B_vec);
iter = 0;
%fprintf('Iteration %d: relative residual = %e\n', iter, res_norm);

%--- Main Conjugate Gradient loop ---
while (res_norm > tol) && (iter < maxIter)
    Ap = applyOperator(p);
    alpha = (r' * r) / (p' * Ap);
    
    % Update solution
    J_vec = J_vec + alpha * p;
    
    % Update residual
    r_new = r - alpha * Ap;
    
    % Compute beta for the next direction
    beta = (r_new' * r_new) / (r' * r);
    
    % Update direction vector
    p = r_new + beta * p;
    
    % Update residual and iteration counter
    r = r_new;
    iter = iter + 1;
    res_norm = norm(r) / norm(B_vec);
    %fprintf('Iteration %d: relative residual = %e\n', iter, res_norm);
end

toc;

% matlab CG
tic;

x = cgs(applyOperator, B_vec, tol, maxIter);

toc;

tic;

x_exact = A \ B_vec;

toc;

norm(J_vec)
norm(x)
norm(x_exact)
