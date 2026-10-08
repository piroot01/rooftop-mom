% optimal_omega_circulant_fixed.m
% Corrected version: builds Chan circulant and optimal omega-circulant preconditioners
% for the tridiagonal Toeplitz T with t0 = 2+alpha, t_{±1} = -1.
% Compares spectra of C^{-1}T and runs cgs, reporting iteration counts.

clear; close all; rng(1);

%% parameters
n = 8;
tol = 1e-6;
maxit = 1000;
b = ones(n, 1);

%% build Toeplitz tridiagonal T
t = zeros(n, 1);
t(1) = 2;
t(2) = -1;

T = toeplitz(t, t);

%% Eigenvalues of the matrix
Tlambda = 2 - 2 * cos(pi * (1:n).' ./ (n + 1));

% first column and first row (as vectors)
col = full(T(:,1));    % col(k+1) = t_k for k >= 0
row = full(T(1,:)).';  % row(k+1) = t_{-k} for k >= 0

%% Chan circulant (uses safe get_t)
C_ch = chan_circulant_from_toeplitz(n, col, row);   % omega=1 inside

%% Compute optimal omega parameter
% s = sum_{j=1}^{n-1} conj(t_{-j}) * t_{n-j}
s = 0;
for j = 1:n-1
    tj_neg = get_t_from_colrow(col, row, -j);
    t_n   = get_t_from_colrow(col, row, n-j);
    s = s + conj(tj_neg) * t_n;
end

if abs(s) < 1e-16
    phi = 0;
else
    phi = angle(s);
end
omega = exp(1i * phi / n)    % omega = exp(i * phi / n)

% build diagonal D
omega_pows = omega .^ (0:n-1);
D = spdiags(omega_pows.', 0, n, n);

% circulant for modified Toeplitz T(omega)
Cn_omega = chan_circulant_from_toeplitz(n, col, row, omega);

% assemble omega-circulant
C_omega = D * Cn_omega * conj(D);

%% Convert to full for eigen computations (moderate n)
C_ch = full(C_ch);
C_omega = full(C_omega);
Tfull = full(T);

%% Regularize tiny nonpositive eigenvalues (if any)
% A very small bump to avoid singulars in backslash, if needed.
eps_reg = 1e-14;
minCch = min(real(eig(C_ch)));
minCow = min(real(eig(C_omega)));
if minCch <= eps_reg
    C_ch = C_ch + (abs(minCch) + eps_reg) * eye(n);
end
if minCow <= eps_reg
    C_omega = C_omega + (abs(minCow) + eps_reg) * eye(n);
end

%% Spectral comparisons
A_ch = C_ch \ Tfull;
A_w  = C_omega \ Tfull;

ev_ch = eig(A_ch);
ev_w  = eig(A_w);

figure(1); clf;
subplot(1,2,1)
plot(real(ev_ch), imag(ev_ch), 'o', 'MarkerSize',6, 'MarkerFaceColor','b');
xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title('Spectrum of C_{Chan}^{-1} T'); grid on; axis equal;

subplot(1,2,2)
plot(real(ev_w), imag(ev_w), 'o', 'MarkerSize',6, 'MarkerFaceColor','r');
xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title('Spectrum of C_{\omega}^{-1} T'); grid on; axis equal;

fprintf('omega parameter: phi = %.12g, omega = %.12g + %.12gi\n', phi, real(omega), imag(omega));
fprintf('Chan eigen real-range: [%g, %g]\n', min(real(ev_ch)), max(real(ev_ch)));
fprintf('omega eigen real-range: [%g, %g]\n', min(real(ev_w)), max(real(ev_w)));

%% Run CGS with left preconditioners
fprintf('\nRunning cgs (tol=%g, maxit=%d)\n', tol, maxit);

[x0, f0, rel0, iter0, resvec0] = cgs(Tfull, b, tol, maxit);
fprintf('Unpreconditioned cgs: flag=%d, iter=%d\n', f0, iter0);

[x_ch, f_ch, rel_ch, iter_ch, resvec_ch] = cgs(Tfull, b, tol, maxit, C_ch, []);
fprintf('Chan preconditioned cgs: flag=%d, iter=%d\n', f_ch, iter_ch);

[x_w, f_w, rel_w, iter_w, resvec_w] = cgs(Tfull, b, tol, maxit, C_omega, []);
fprintf('omega-circ preconditioned cgs: flag=%d, iter=%d\n', f_w, iter_w);

%% Plot residual histories
figure(2); clf;
semilogy(0:length(resvec0)-1, resvec0 / norm(b), '-d','LineWidth',1.5); hold on;
semilogy(0:length(resvec_ch)-1, resvec_ch / norm(b), '-x','LineWidth',1.5);
semilogy(0:length(resvec_w)-1,  resvec_w / norm(b),  '-o','LineWidth',1.5);
legend('unpreconditioned','Chan','\omega-circ','Location','best');
xlabel('Iteration'); ylabel('Relative residual ||r_k|| / ||b||');
title('CGS residual history');
grid on;

fprintf('\nSummary (n=%d, alpha=%g):\n', n, alpha);
fprintf('  Chan:   real-eig range [%.4e, %.4e], cgs iters = %d\n', min(real(ev_ch)), max(real(ev_ch)), iter_ch);
fprintf('  omega:  real-eig range [%.4e, %.4e], cgs iters = %d\n', min(real(ev_w)), max(real(ev_w)), iter_w);


%% -----------------------------
%% Local helper functions below
%% -----------------------------

function C = chan_circulant_from_toeplitz(n, col, row, omega)
% Build Chan circulant using t_j (from col,row). Optional omega applies
% t_j -> t_j * omega^j.
if nargin < 4
    omega = 1;
end
c = zeros(n,1);
for k = 0:n-1
    % t_k and t_{k-n} (safe access via get_t_from_colrow)
    tk    = get_t_from_colrow(col, row, k) * (omega^k);
    tk_n  = get_t_from_colrow(col, row, k-n) * (omega^(k-n));
    c(k+1) = ( (n-k)*tk + k*tk_n ) / n;
end
% turn c (length-n vector c_k, k=0..n-1) into first column of circulant
first_col = [c(1); flipud(c(2:end))];
first_row = c.';
C = toeplitz(first_col, first_row);
end

function val = get_t_from_colrow(col, row, j)
% Safe retrieval of Toeplitz coefficient t_j using first column and first row.
% col(k+1) = t_k for k >= 0 (0..n-1)
% row(k+1) = t_{-k} for k >= 0 (0..n-1)
ncol = numel(col);
nrow = numel(row);
if j >= 0
    idx = j + 1;
    if idx >= 1 && idx <= ncol
        val = col(idx);
    else
        val = 0;
    end
else
    idx = -j + 1;
    if idx >= 1 && idx <= nrow
        val = row(idx);
    else
        val = 0;
    end
end
end
