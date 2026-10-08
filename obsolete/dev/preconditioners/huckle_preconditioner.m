function demo_toeplitz_preconditioners()
    % DEMO_TOEPLITZ_PRECONDITIONERS
    %
    % This script demonstrates and compares circulant preconditioners 
    % for Hermitian Toeplitz matrices, including:
    % 1. T. Chan's Optimal Preconditioner (CF)
    % 2. Tyrtyshnikov's Superoptimal Preconditioner (CT) - Using Fast O(N log N) construction
    % 3. Huckle's Frobenius Preconditioner (CB)
    %
    % It compares these against the CGS solver and standard PCG.

    close all; clc; clear;

    % --- Configuration ---
    N = 64; 
    tol = 1e-8;
    max_it = 100;
    
    fprintf('Running Tests for N = %d\n', N);
    fprintf('------------------------------------------------------------\n');

    % --- Define Test Cases from Tyrtyshnikov (1992) ---
    test_cases = {};
    
    % Tyrtyshnikov Ex 1: Tridiagonal [2, -1]
    col1 = zeros(N,1); col1(1)=2; col1(2)=-1;
    test_cases{end+1} = struct('name', 'Tyrtyshnikov Ex1 (Tridiag)', 'col', col1);
    
    % Tyrtyshnikov Ex 2: 1 / (0.25 + (i-j)^2)
    idx = (0:N-1)';
    col2 = 1 ./ (0.25 + idx.^2);
    test_cases{end+1} = struct('name', 'Tyrtyshnikov Ex2 (Rational)', 'col', col2);
    
    % Tyrtyshnikov Ex 3: 1 / (1 + |i-j|)
    col3 = 1 ./ (1 + idx);
    test_cases{end+1} = struct('name', 'Tyrtyshnikov Ex3 (Slow Decay)', 'col', col3);

    % Huckle Example (c)
    col4 = zeros(N,1); col4(1)=0.6;
    idx_h = (1:N-1)';
    col4(2:end) = cos(idx_h) ./ (idx_h + 1);
    test_cases{end+1} = struct('name', 'Huckle Ex (c)', 'col', col4);

    % --- Main Loop ---
    for t = 1:length(test_cases)
        tc = test_cases{t};
        fprintf('\nTest Case: %s\n', tc.name);
        
        c = tc.col; % Symmetric Toeplitz defined by first column
        
        % RHS
        b = ones(N, 1);
        b = b / norm(b);

        % --- 1. Construct Preconditioners ---
        
        % A) T. Chan's Optimal (CF)
        tic;
        cf_eig = build_tchan_fast(c);
        time_cf = toc;
        
        % B) Tyrtyshnikov's Superoptimal (CT)
        % Uses purely FFT based construction (O(N log N)) without building T
        tic;
        ct_eig = build_tyrtyshnikov_fast(c);
        time_ct = toc;
        
        % C) Huckle's Frobenius (CB)
        % Uses FFTs to build the system D*x=b
        tic;
        cb_eig = build_huckle_fast(c);
        time_cb = toc;

        % --- 2. Iterative Solution ---
        % Matrix-vector product function using FFT
        afun = @(x) fast_toeplitz_sym_mult(c, x);
        
        % PCG with Preconditioners
        [~,fl0,rr0,it0,rv0] = pcg(afun, b, tol, max_it);
        [~,fl1,rr1,it1,rv1] = pcg(afun, b, tol, max_it, @(x) solve_circ(cf_eig, x));
        [~,fl2,rr2,it2,rv2] = pcg(afun, b, tol, max_it, @(x) solve_circ(ct_eig, x));
        [~,fl3,rr3,it3,rv3] = pcg(afun, b, tol, max_it, @(x) solve_circ(cb_eig, x));
        
        % CGS (Conjugate Gradient Squared) - No Preconditioner
        % Note: CGS is for non-symmetric usually, but requested for comparison
        [~,fl_cgs,rr_cgs,it_cgs,rv_cgs] = cgs(afun, b, tol, max_it);

        % --- Report ---
        fprintf('%-15s | %-8s | %-5s | %-9s\n', 'Solver', 'Time(s)', 'Iters', 'RelRes');
        fprintf('--------------------------------------------------\n');
        fprintf('%-15s | %-8.4f | %-5d | %-9.2e\n', 'PCG (None)', 0, it0, rr0);
        fprintf('%-15s | %-8.4f | %-5d | %-9.2e\n', 'PCG (T.Chan)', time_cf, it1, rr1);
        fprintf('%-15s | %-8.4f | %-5d | %-9.2e\n', 'PCG (Tyrt)', time_ct, it2, rr2);
        fprintf('%-15s | %-8.4f | %-5d | %-9.2e\n', 'PCG (Huckle)', time_cb, it3, rr3);
        fprintf('%-15s | %-8s | %-5d | %-9.2e\n', 'CGS (None)', '-', length(rv_cgs)-1, rr_cgs);

        % --- Spectra Analysis ---
        % Reconstruct full matrix just for eig() calculation to plot spectra
        T_full = toeplitz(c);
        
        % Eigenvalues of C^-1 T
        lam_T = eig(T_full);
        lam_CF = eig(T_full * inv(circ_from_eig(cf_eig))); 
        lam_CT = eig(T_full * inv(circ_from_eig(ct_eig)));
        lam_CB = eig(T_full * inv(circ_from_eig(cb_eig)));
        
        figure('Name', ['Spectra: ' tc.name]);
        subplot(2,2,1); plot(sort(real(lam_T)), 'k.'); title('Original T'); grid on;
        subplot(2,2,2); plot(sort(real(lam_CF)), 'b.'); title('T. Chan C_F^{-1}T'); grid on;
        subplot(2,2,3); plot(sort(real(lam_CT)), 'r.'); title('Tyrtyshnikov C_T^{-1}T'); grid on;
        subplot(2,2,4); plot(sort(real(lam_CB)), 'g.'); title('Huckle C_B^{-1}T'); grid on;
    end
end

%% --- Fast Constructors using FFTs ---

function eig_vals = build_tchan_fast(c)
    % T. Chan optimal preconditioner from first column c
    % O(N) construction
    n = length(c);
    c_chan = zeros(n, 1);
    c_chan(1) = c(1);
    k = (1:n-1)';
    % Formula: c_k = ((n-k)t_k + k t_{n-k}) / n
    % Since symmetric t_k = t_{n-k} ? No, t is first column.
    % T is symmetric => row is c'. 
    % diag k of T has value c(k+1). diag -(n-k) has value c(n-k+1).
    % Formula for symmetric T:
    c_chan(2:end) = ((n-k).*c(2:end) + k.*c(end:-1:2)) / n;
    eig_vals = fft(c_chan);
end

function eig_vals = build_tyrtyshnikov_fast(c)
    % Tyrtyshnikov Superoptimal (C_T)
    % Solves (U + U^T) c_inv = f
    % U = circulant matrix of traces of T*Q^k*T
    % Uses O(N log N) method from Section 5 (Theorem 5.1)
    
    n = length(c);
    
    % 1. Compute f vector (traces of cyclic diagonals)
    % f_k = 2 * sum(k-th cyclic diag of T)
    % For symmetric T, cyclic diag k has elements c(k+1) and c(n-k+1)
    k = (0:n-1)';
    % diag k has n-k elements of c(k+1) and k elements of c(n-k+1)
    % logic: index 0 is main diag.
    c_rev = [c(1); c(end:-1:2)];
    f = 2 * ((n-k).*c + k.*c_rev);
    
    % 2. Compute first column of U (vector of traces s)
    % u_k = tr(T * Q^k * T). 
    % We compute s(T^2) where s_k(M) is sum of k-th cyclic diag of M.
    % Decomposition T = L + U_mat. 
    % Let L_vec = c (first col). L matrix is lower triangular Toeplitz.
    % Let R_vec = [0; c(2:end)]. R matrix is strict upper triangular Toeplitz.
    % T = L + R.
    % T^2 = L^2 + LR + RL + R^2.
    
    % Term 1: L^2 (Lower * Lower)
    % First col of L^2 is conv(c, c) truncated
    l2_col = ifft(fft(c, 2*n) .^ 2);
    l2_col = l2_col(1:n);
    % Sum of cyclic diags of L^2. 
    % L^2 is lower triangular. k-th cyclic diag has (n-k) elements from l2_col(k+1)
    % and k elements from 0 (upper part is 0).
    s_LL = (n - k) .* l2_col;
    
    % Term 4: R^2 (Strict Upper * Strict Upper)
    % Similar to L^2 but shifted. R is L without diagonal.
    r_vec = [0; c(2:end)];
    r2_row = ifft(fft(r_vec, 2*n) .^ 2); % This gives first col of "Lower" version
    r2_row = r2_row(1:n);
    % R^2 is upper. k-th cyclic diag has k elements of value?
    % Mapping: Upper diag k corresponds to cyclic index n-k (or -k).
    % Let's use symmetry. s(R^2) should be reverse of s(L_strict^2).
    % Actually simpler: U is symmetric if T is symmetric.
    % So s(T^2) is symmetric vector.
    % We only need to compute s(L^2 + LR). The rest follows by symmetry?
    % T^2 is symmetric => s_k(T^2) = s_{n-k}(T^2).
    
    % Term 2: LR (Lower * Strict Upper)
    % Use Tyrtyshnikov Theorem 5.1
    % M = L_mat * R_mat.
    % l_vec = c. r_vec_row = [0, c(2), ..., c(n)].
    r_vec_row = [0; c(2:end)];
    s_LR = fast_trace_LR(c, r_vec_row);
    
    % Term 3: RL (Strict Upper * Lower)
    % Use Corollary 2: s(RL) can be computed via J*RL*J
    % J*R*L*J = (J*R*J) * (J*L*J).
    % J*R*J is Strict Lower. J*L*J is Upper.
    % This is again form Lower * Upper.
    % l_new = J * r_vec_row (as col). r_new = J * c (as row).
    l_new = flipud(r_vec_row); % Strict Lower
    r_new = flipud(c);         % Upper
    % Wait, J*T*J flips matrix.
    % J*R*J is a Toeplitz matrix? Yes.
    % If R is Toeplitz(0, [0, t1...]), J*R*J is Toeplitz([0, t_n-1...], 0).
    % Effectively, we compute s(Lower_new * Upper_new).
    s_RL_raw = fast_trace_LR(l_new, r_new);
    % s_k(M) = s_{n-k}(J M J).
    s_RL = [s_RL_raw(1); s_RL_raw(end:-1:2)];
    
    % Total s vector
    s_total = s_LL + s_LR + s_RL + [0; s_LL(end:-1:2)]; % R^2 is symmetric to L^2? 
    % Correction for R^2:
    % R corresponds to L with 0 diagonal.
    % R^2 is strict upper.
    % s_k(R^2) involves upper diags. Upper diag k is cyclic n-k.
    % s(R^2) should be reflection of s( (L-diag)^2 ).
    l_strict = c; l_strict(1) = 0;
    l2_strict_col = ifft(fft(l_strict, 2*n).^2); l2_strict_col = l2_strict_col(1:n);
    s_LL_strict = (n-k).*l2_strict_col;
    s_RR = [s_LL_strict(1); s_LL_strict(end:-1:2)]; % Flip for upper
    
    s_total = s_LL + s_LR + s_RL + s_RR;
    
    % This s_total is the first column of U
    
    % 3. Solve (2*U) * hat_c = f
    % U is circulant with first col s_total
    lam_U = fft(s_total);
    lam_f = fft(f);
    
    % Check for singularity (rare for positive definite T)
    lam_inv_c = (lam_f / 2) ./ lam_U;
    
    % Result is eigenvalues of C^-1.
    % We want eigenvalues of C.
    eig_vals = 1 ./ lam_inv_c;
    
    % Enforce real
    eig_vals = real(eig_vals);
end

function s = fast_trace_LR(l, r)
    % Implements Eq 5.7 from Tyrtyshnikov (1992)
    % Computes traces s_k of product L * R where L is lower Toeplitz (col l)
    % and R is upper Toeplitz (row r).
    n = length(l);
    
    % Vectors for P and Q matrix-vector mults
    % P and Q are Toeplitz. We compute P*v and Q*w using FFT padding.
    
    % Prepare vectors v1, v2 for P and Q
    % From Eq 5.7:
    % Input to P: [r_n-1 ... r_0]'. In our notation r is [r_0 ... r_n-1].
    % So v_P = flipud(r).
    v_P = flipud(r(:));
    
    % Input to Q: [0, r_n-1, 2r_n-2, ..., (n-1)r_1]'
    % indices of r: 1..n corresponds to r0..rn-1
    % multiplier array: 0..n-1
    mults = (0:n-1)';
    v_Q = mults .* flipud(r(:)); % Check indexing carefully
    % Eq 5.7 vector: [0; r_{n-1}; 2 r_{n-2}; ... ]
    % My v_Q: 0*r_{n-1} ... NO.
    % Matrix Q definition (5.8):
    % Cols are scaled. 
    % Let's use the exact vector definition from (5.7) RHS:
    % Vector 2: [0, r_{n-1}, 2r_{n-2}, ..., (n-1)r_1]^T.
    vec2 = zeros(n,1);
    for k=1:n-1
       vec2(k+1) = k * r(n-k+1); % r(index) is r_{index-1}
       % k=1: 1 * r_{n-1}. Correct.
    end
    
    % P matrix (5.8):
    % Toeplitz. First col: [l_{n-1} ... l_0] * scale? No.
    % P = [ l_{n-1}  2l_{n-2} ...
    %       0        l_{n-1} ... ]
    % P is Upper Triangular Toeplitz?
    % Look at (5.8). P is Upper Toeplitz.
    % First row: l_{n-1}, 2l_{n-2}, ..., n l_0.
    % First col: l_{n-1}, 0, ..., 0.
    % Definition of P: P_{ij} = (j-i+1) l_{n-1-(j-i)}. (for j>=i)
    
    % Q matrix (5.8):
    % Lower Triangular Toeplitz.
    % First col: l_0, l_1, ..., l_{n-1}. This is just L!
    % Wait, Q = [ l_0   0
    %             l_1   l_0
    %             ...       ]
    % Yes, Q is L.
    % Wait, Q definition in 5.8:
    % Q = [ l_0 0 ...
    %       l_1 l_0 ... ]
    % So Q is exactly the matrix L defined by vector l.
    
    % Compute P * v_P
    % P is Upper Toeplitz defined by row p_row.
    % p_row = [ l_{n-1}, 2 l_{n-2}, ..., n l_0 ].
    p_row = zeros(1, n);
    for k=0:n-1
        p_row(k+1) = (k+1) * l(n-k); % l(n-k) is l_{n-1-k}
    end
    % Mult Upper Toeplitz P with v_P
    % Embed in circulant 2n.
    % Col for circulant: [P_{00}, 0...0, P_{0n}...P_{01}]' -> [p_row(1), 0..0, p_row(n)..p_row(2)]
    p_circ_col = [p_row(1); zeros(n,1); p_row(end:-1:2)'];
    v_P_pad = [v_P; zeros(n,1)];
    res_P = ifft(fft(p_circ_col) .* fft(v_P_pad));
    res_P = res_P(1:n);
    
    % Compute Q * vec2
    % Q is Lower Toeplitz defined by col l.
    % Col for circulant: [l; 0...; 0] (standard Toeplitz embed)
    q_circ_col = [l; zeros(n,1)]; % l is [l0...ln-1]
    vec2_pad = [vec2; zeros(n,1)];
    res_Q = ifft(fft(q_circ_col) .* fft(vec2_pad));
    res_Q = res_Q(1:n);
    
    s = res_P + res_Q;
end

function eig_vals = build_huckle_fast(c)
    % Huckle C_B preconditioner
    % Solves D * x = b, where D = B .* conj(B)
    % B = F * T * F'
    % x = 1./eig_vals
    
    n = length(c);
    
    % 1. Construct B column by column using FFTs (O(N^2 log N))
    % Col j of B = F * (T * (F' * e_j))
    % T is symmetric Toeplitz.
    B = zeros(n, n);
    I = eye(n);
    
    % Precompute FFT of T's embedding for fast multiply
    c_embed = [c; 0; c(end:-1:2)];
    fft_T = fft(c_embed);
    
    for j = 1:n
        % ej in Fourier domain (column of F')
        % F' * e_j is just the j-th column of F' (inverse FFT matrix * n?)
        % F_{km} = w^(-km). F'_{km} = w^(km).
        % Let's use ifft for F'
        v = I(:, j);
        fv = ifft(v) * n; % F' * v? 
        % Matlab F is not normalized? 
        % definition: B = F T F'.
        % If we use normalized F (unitary), B is unitary transform.
        % Let's assume Unitary DFT for the preconditioner definition.
        % F_norm = fft(I) / sqrt(n).
        % B = F_norm * T * F_norm'.
        
        % Col j of F_norm': 
        % u_j = F_norm' * e_j. 
        % This is inverse FFT of e_j * sqrt(n).
        u_j = ifft(v) * sqrt(n);
        
        % T * u_j using fast mult
        u_pad = [u_j; zeros(n,1)];
        Tu_j = ifft(fft_T .* fft(u_pad));
        Tu_j = Tu_j(1:n);
        
        % F_norm * (T * u_j)
        % This is fft(Tu_j) / sqrt(n)
        col_B = fft(Tu_j) / sqrt(n);
        
        B(:, j) = col_B;
    end
    
    % 2. Construct D = B .* conj(B)
    D = B .* conj(B);
    
    % 3. RHS b = diag(B)
    rhs = diag(B);
    
    % 4. Solve D x = rhs
    % D is positive definite
    x = D \ rhs;
    x = abs(x); % Safety
    
    eig_vals = 1 ./ x;
end

%% --- Utilities ---

function y = fast_toeplitz_sym_mult(c, x)
    % Fast multiplication T*x for symmetric Toeplitz T defined by col c
    n = length(c);
    v = [c; 0; c(end:-1:2)];
    x_pad = [x; zeros(n+1, 1)];
    res = ifft(fft(v) .* fft(x_pad));
    y = res(1:n);
end

function y = solve_circ(eig_vals, x)
    % Apply C^-1 x
    y = ifft(fft(x) ./ eig_vals);
end

function C = circ_from_eig(eig_vals)
    % Reconstruct circulant matrix from eigenvalues
    c = ifft(eig_vals);
    C = toeplitz(c, [c(1) c(end:-1:2).']);
end