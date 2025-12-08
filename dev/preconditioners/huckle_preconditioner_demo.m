clc; clearvars; close all;

%% Build the first row of the Toeplitz matrix T
N = 9;

t = zeros(N, 1);
t(1) = 2;
t(2) = -1;

T = toeplitz(t, t);

Tlambda = 2 - 2 * cos(pi * (1:N).' ./ (N + 1))

b = randn(N, 1);

%% Fast Toeplitz matrix applicator
function data = prepFFTMultData(t)
    data.N = numel(t);
    data.cf = fft([t; 0; flipud(t(2:end))]);
end

function y = FFTMult(data, x)
    N = data.N;
    xp = [x; zeros(N, 1)];
    yp = ifft(data.cf .* fft(xp));
    y = yp(1:N);
end

%% T. Chan's preconditioner
% Build the first column c of the circulant.
function c = buildTChanPrecond(t)
    N = numel(t);
    c = zeros(N,1);
    for j = 0:(N - 1)
        if j == 0
            c(1) = t(1);
        else
            c(j+1) = ((N - j) * t(j+1) + j * t(N - j + 1)) / N;
        end
    end
end

%% Optimal {omega}-circulant preconditioner
function val = get_t(t, j)
    N = numel(t);
    if j >= 0
        idx = j + 1;
        if idx >=1 && idx <= N
            val = t(idx);
        else
            val = 0;
        end
    else
        idx = -j + 1;
        if idx >=1 && idx <= N
            val = t(idx);
        else
            val = 0;
        end
    end
end

function omega = findOptimalOmega(t)
    N = numel(t);
    s = 0;
    for n = 1:(N - 1)
        tk = get_t(t, n);
        tk_minus_N = get_t(t, n - N);
        s = s + conj(tk) * tk_minus_N;
    end
    omega = exp(1i * angle(s) / N);
end

omegaOpt = findOptimalOmega(t);

function c = buildOmegaPrecond(t, omega)
    N = numel(t);
    c = zeros(N, 1);
    for n = 0:(N-1)
        tk_omega = get_t(t, n) * omega^n;
        tk_minus_N_omega = get_t(t, n - N) * omega^(n - N);
        % Chan-style entry
        c(n + 1) = ((N - n) * tk_omega + n * tk_minus_N_omega) / N;
    end
end

%% Circulant preconditioner applicator using FFT
function y = FFTPrecond(pf, x)
    y = ifft(fft(x) ./ pf);
end

%% Direct spectra knowledge based preconditioner
% --- DST-I using FFT (returns S * x where S_{j,k} = sin(pi*j*k/(n+1)) ) ---
function y = dst1(x)
    n = numel(x);
    % build odd extension and use FFT
    X = fft([0; x(:); 0; -flipud(x(:))]);
    % indices: MATLAB 1-based, take elements 2..n+1, multiply by -1/2
    y = -0.5 * imag(X(2:n+1));
end

% inverse of the (unnormalized) DST-I is the same DST-I up to sign/scale,
% but we will handle the orthonormalization below.

% --- apply T^{-1} via sine-basis (fast) ---
function y = apply_Tinv_via_DST(b, lambda)
    % lambda : vector of eigenvalues (length n)
    n = numel(b);
    % normalization factor for orthonormal sine basis: sqrt(2/(n+1))
    alpha = sqrt(2/(n+1));
    % forward transform (S' * b). For our symmetry DST-I gives S'*b = S*b
    scoeff = alpha * dst1(b);        % = S_mat' * b
    % diagonal solve in spectral domain
    scoeff = scoeff ./ lambda(:);
    % inverse transform
    y = alpha * dst1(scoeff);        % = S_mat * (Lambda^{-1} * S_mat' * b)
end

%% Prepare the data for FFT multiplication and construct the handle
dataFFTMult = prepFFTMultData(t);
dataFFTMult.cf
FFTMultHandle = @(x) FFTMult(dataFFTMult, x);

%% Direct solver

xRef = T \ b;

%% Conjugate gradient solver
tolerance = 1e-12;
maxIterCount = 1000;

% Unpreconditioned variant - spatial domain
t0 = tic;
[x0, flag0, relres0, iter0, resvec0] = cgs(FFTMultHandle, b, tolerance, maxIterCount);
solverTime = toc(t0);
fprintf('Unpreconditioned: flag = %d, iteration count = %d, solver time = %g s.\n', flag0, iter0, solverTime);

% Preconditioned variant - T. Chan's circulant preconditioner
tChanPrecond = buildTChanPrecond(t);
tChanLambda = fft(tChanPrecond);
applyHandle = @(x) FFTPrecond(tChanLambda, FFTMultHandle(x));

t0 = tic;
[x_TChan, flag1, relres1, iter1, resvec_TChan] = cgs(applyHandle, FFTPrecond(tChanLambda, b), tolerance, maxIterCount);
solverTime = toc(t0);
fprintf('Unpreconditioned: flag = %d, iteration count = %d, solver time = %g s.\n', flag1, iter1, solverTime);

% Preconditioned variant - optimal {omega}-circulant
omegaPrecond = buildOmegaPrecond(t, omegaOpt);
omegaLambda = fft(omegaPrecond);
applyHandle = @(x) FFTPrecond(omegaLambda, x);

t0 = tic;
[x_omega, flag2, relres2, iter2, resvec_Omega] = cgs(FFTMultHandle, b, tolerance, maxIterCount, applyHandle);
solverTime = toc(t0);
fprintf('Unpreconditioned: flag = %d, iteration count = %d, solver time = %g s.\n', flag2, iter2, solverTime);

% Preconditioned variant - known spectra of T
applyHandle = @(x) apply_Tinv_via_DST(x, Tlambda);

t0 = tic;
[xTlambda, flag3, relres3, iter3, resvec_Tlambda] = cgs(@(x) applyHandle(FFTMultHandle(x)), applyHandle(b), tolerance, maxIterCount);
solverTime = toc(t0);
fprintf('Unpreconditioned: flag = %d, iteration count = %d, solver time = %g s.\n', flag3, iter3, solverTime);

norm(xRef - xTlambda)

%% Residual vector convergence comparison
figure;
clf;
semilogy(0:length(resvec0)-1, resvec0 / norm(b), '-o','LineWidth',1.5);
hold on;
semilogy(0:length(resvec_TChan)-1, resvec_TChan / norm(b), '-x','LineWidth',1.5);
semilogy(0:length(resvec_Omega)-1,  resvec_Omega / norm(b),  '-d','LineWidth',1.5);
semilogy(0:length(resvec_Tlambda)-1,  resvec_Tlambda / norm(b),  '-d','LineWidth',1.5);
legend('unpreconditioned','Chan','\omega-circ', 'T lambda directly','Location','best');
xlabel('Iteration'); ylabel('Relative residual ||r_k|| / ||b||');
title('CGS residual history');
grid on;