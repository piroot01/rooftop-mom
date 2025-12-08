% demo_toeplitz_cg_minfft.m
% Demo of PCG for symmetric Toeplitz using minimal FFT strategy
% (based on Huckle, "Iterative Methods for Toeplitz-like Matrices",
%  representation (8) and discussion of reduced FFTs). :contentReference[oaicite:1]{index=1}
clear; close all; clc;

%% Problem setup: symmetric Toeplitz (e.g. discretized 1D Laplacian + small diag)
n = 1024;                  % size (choose power of 2 for speed)
kappa = 1e-2;              % small regularization
j = (0:n-1).';
% construct first column of symmetric Toeplitz, example: Laplacian-like kernel
t = zeros(n,1);
t(1) = 2 + kappa;
t(2) = -1;
% for larger bandwidth use t(3:n) = ... ; here simple tridiagonal Toeplitz
% the Toeplitz matrix is symmetric positive definite for this choice
% Right hand side
x_true = sin((1:n)'*pi/ (n+1));        % a known solution
% compute b = T * x_true using a fast toeplitz multiply routine below
b = toeplitz_matvec_fft(t, x_true);    % uses the same circulant-embedding routine

%% Solve with our PCG that reuses FFTs and uses a circulant preconditioner
tol = 1e-10;
maxit = 200;
fprintf('Running PCG (FFT-accelerated Toeplitz matvec + circulant preconditioner)...\n');
tic
[x_pcg, flag, relres, iter, resvec] = toeplitz_pcg_minfft(t, b, tol, maxit);
t_pcg = toc;
fprintf('PCG done: flag=%d, iter=%d, relres=%.2e, time=%.3fs\n', flag, iter, relres, t_pcg);

%% Compare with MATLAB pcg on explicit Toeplitz (only for verification / small n)
if n <= 2000
    T = toeplitz(t);
    tic
    [x_matlab, f2, r2, it2] = pcg(T, b, tol, maxit);
    t_matlab = toc;
    fprintf('MATLAB pcg on explicit Toeplitz: flag=%d, iter=%d, time=%.3fs\n', f2, it2, t_matlab);
    err = norm(x_pcg - x_matlab)/norm(x_matlab);
    fprintf('Relative difference between solvers: %.2e\n', err);
end

%% Plot residual history
figure;
semilogy(0:iter, resvec / resvec(1), '-o');
xlabel('Iteration'); ylabel('relative residual');
title('PCG residual history (normalized)');

%% Report solution error against true x
err_true = norm(x_pcg - x_true)/norm(x_true);
fprintf('Relative error vs true solution: %.3e\n', err_true);

%% ---------------- Local functions ----------------

function y = toeplitz_matvec_fft(t, x)
% Toeplitz matrix-vector product using circulant embedding and FFT.
% Input:
%   t : first column of symmetric Toeplitz (n x 1)
%   x : vector (n x 1)
% Output:
%   y = T * x (n x 1)
    n = numel(t);
    % choose m = power of two >= 2*n-1 for fastest FFT
    m = 2^nextpow2(2*n-1);
    % build first column of circulant embedding size m
    tcirc = zeros(m,1);
    tcirc(1:n) = t;
    % place t(2:n) reversed at tail positions to complete circulant first column
    tcirc(m - (1:(n-1)) + 1) = t(2:n);
    % convolution via FFT
    X = fft([x; zeros(m-n,1)]);
    Y = fft(tcirc) .* X;
    ybig = ifft(Y);
    y = real(ybig(1:n));   % numerical rounding may produce tiny imag
end

function [x,flag,relres,iter,resvec] = toeplitz_pcg_minfft(t,b,tol,maxit)
% PCG solver for symmetric Toeplitz T defined by first column t.
% Implements matvec by circulant embedding FFT and uses an (optimal-like)
% circulant preconditioner solved via FFT diagonalization (all precomputed).
%
% Outputs follow MATLAB's pcg style:
%   x : approximate solution
%   flag : 0 if converged, 1 if not converged after maxit
%   relres : final relative residual
%   iter : # iterations used
%   resvec : residual norm history (length iter+1)
%
    n = numel(t);
    m = 2^nextpow2(2*n-1);         % embed length (power of 2 recommended)

    % ---------- Build circulant embedding of T (for matvec) ----------
    tcirc = zeros(m,1);
    tcirc(1:n) = t;
    tcirc(m - (1:(n-1)) + 1) = t(2:n);
    Tfft = fft(tcirc);             % precompute FFT of circulant kernel

    % ---------- Build circulant preconditioner (Chan-like for symmetric) ----------
    % For symmetric Toeplitz with first column t0,t1,... we can set
    % t_neg(k) ~ t_{k-n} which for symmetric == t(n-k)
    t_neg = [ t(1); flipud(t(2:end)) ];    % length n, approx. t_{k-n}
    c = zeros(n,1);
    for k = 0:(n-1)
        c(k+1) = ((n-k)/n) * t(k+1) + (k/n) * t_neg(k+1);
    end
    % Build circulant embedding for preconditioner (length m) and its FFT
    ccirc = zeros(m,1);
    ccirc(1:n) = c;
    ccirc(m - (1:(n-1)) + 1) = c(2:n);
    Cfft = fft(ccirc);

    % -- safeguard small entries in Cfft (avoid division by zero) --
    eps_reg = 1e-12 * max(abs(Cfft));
    Cfft(abs(Cfft) < eps_reg) = eps_reg;

    % ---------- helper nested functions that use precomputed FFTs ----------
    function y = Aprod(xv)
        % compute y = T * xv using circulant embedding (size m)
        X = fft([xv; zeros(m-n,1)]);
        Y = Tfft .* X;
        ybig = ifft(Y);
        y = real(ybig(1:n));
    end

    function z = Minv_apply(r)    % apply circulant preconditioner (solve C z = r)
        % Solve C z = r using diagonalization of circulant ccirc (FFT)
        R = fft([r; zeros(m-n,1)]);
        Z = R ./ Cfft;
        zbig = ifft(Z);
        z = real(zbig(1:n));
    end

    % ---------- PCG initialization ----------
    x = zeros(n,1);
    r = b - Aprod(x);
    normb = norm(b);
    if normb == 0
        normb = 1;
    end
    res0 = norm(r);
    resvec = res0;
    if res0 <= tol * normb
        flag = 0; relres = res0 / normb; iter = 0; return;
    end

    z = Minv_apply(r);
    p = z;
    rz_old = r' * z;

    flag = 1;
    iter = 0;
    for k = 1:maxit
        Ap = Aprod(p);
        alpha = rz_old / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        res = norm(r);
        resvec(end+1) = res; %#ok<AGROW>
        relres = res / normb;
        iter = k;
        if (res <= tol * normb)
            flag = 0;
            break;
        end
        z = Minv_apply(r);
        rz_new = r' * z;
        beta = rz_new / rz_old;
        p = z + beta * p;
        rz_old = rz_new;
    end
    relres = norm(r) / normb;
end
