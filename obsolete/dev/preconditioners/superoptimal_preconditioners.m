clc;

% Parameters
%{
n = 10; % Matrix size
ht = rand([n, 1])
%}
n = 16;
ht = zeros(n, 1);
ht(1) = 2.007;
ht(2) = -1;

b = ones(n, 1);

T = toeplitz(ht, ht);

c2 = zeros(n, 1);

T2 = T*T;

for k = 1:n
    ck = zeros(n, 1);
    for j = 1:n
        ck(j) = ck(j) + T2(j, mod(j + k, n) + 1);
    end
    c2(k) = 1 / n *sum(ck);
end

lambda_2 = fft(c2);

% --- 1. T. Chan's Optimal Preconditioner (C_F) ---
% Formula: c_k = ((n-k)t_k + k*t_{n-k})/n
c_chan = zeros(n,1);
for j = 0:(n-1)
    if j == 0
        c_chan(1) = ht(1);   % j = 0 => c0 = a0
    else
        % a_j is a(j+1); a_{j-n} = a_{n-j} = a(N-j+1)
        c_chan(j+1) = ((n - j) * ht(j+1) + j * ht(n - j + 1)) / N;
    end
end

lambda_chan = fft(c_chan); % Eigenvalues of C_F

lambda_tyr = lambda_2 ./ lambda_chan;

apply_prec_ = @(x) apply_prec(x, lambda_tyr);

tol = 1e-6;
maxIt = 100;
[x, flag, ~, iter] = cgs(T, b, tol, maxIt, apply_prec_)
eigs(@(x) apply_prec_(T * x), n, n)
%eigs(T, n)

%{
% --- 2. Tyrtyshnikov's Superoptimal Preconditioner (C_T) ---
% Formula: C_T = Circ(T^2) * Circ(T)^(-1)
% Note: While this is the definition, the paper mentions C_F2 = Circ(T^2)^(1/2) 
% is often preferred. However, adhering to the standard definition of C_T:

% Calculate first column of T^2 efficiently or explicitly
T2 = T * T;

col_T2 = T2(:, 1); 

% Compute Circ(T^2) using Chan's formula on the column of T^2
c_T2 = zeros(n, 1);
for j = 0:(n-1)
    if j == 0
        c_T2(1) = col_T2(1);
    else
        c_T2(j+1) = ((n - j) * col_T2(j+1) + j * col_T2(n - j + 1)) / n;
    end
end
lambda_T2 = fft(c_T2);

% Eigenvalues of C_T = lambda(Circ(T^2)) ./ lambda(Circ(T))
lambda_tyr = lambda_T2 ./ lambda_chan;

% --- 3. Huckle's Superoptimal Preconditioner (C_B) ---
% Minimize ||I - C^(-1/2) A C^(-1/2)||_F
% Algorithm requires solving non-linear system: b_kk = sum(|b_kj|^2 / lambda_j)
% The text suggests applying PCG to D*x = b where x = 1./lambda

F = fft(eye(n)) / sqrt(n); % Normalized Fourier matrix
B = F * T * F';            % B = F T F^H
D = B .* conj(B);          % Hadamard product D = B o B_bar
b_diag = diag(B);          % Diagonal of B (which are eigenvalues of C_F roughly)

% We need to solve D * x = b_diag for x, where lambda = 1./x
% Use PCG to solve this smaller system. 
% Preconditioner for this inner PCG is diag(b_11^2, ... b_nn^2) 
M_inner = diag(b_diag.^2);
[x_vec, ~] = pcg(D, b_diag, 1e-6, 20, M_inner);

% Handle potential negative components (Huckle Suggestion)
% "If x has negative components... define lambda_i = |1/x_i|" 
lambda_huckle = abs(1 ./ x_vec);
%}
% Helper to apply circulant preconditioner in O(n log n)
function y = apply_prec(r, eigenvalues)
    % M^{-1} r = F^H * Lambda^{-1} * F * r
    % MATLAB's fft is unnormalized, so we adjust carefully.
    % ifft(fft(r) ./ eigenvalues) is the standard fast apply.
    y = ifft(fft(r) ./ eigenvalues);
end
