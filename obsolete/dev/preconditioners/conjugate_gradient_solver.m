clc; clearvars;

% First we generate random symmetric Toeplitz matrix
%ht = randn([N, 1]);

N = 10;


ht = zeros(N, 1);
for n = 1:N
    ht(n) = (1 + abs(n))^(-1.1);
end

%{
ht = zeros(N, 1);
ht(1) = 2.006;
ht(2) = -1;
%}
%{
t0 = pi^2/3;
tk = @(k) 4 * (-1)^k / k^2; % for k>=1
ht = [t0; arrayfun(@(k) tk(k), 1:N-1).'];
%}

%{
% user parameters (play with these)
N = 200;            % number of basis functions (matrix size NxN)
k = 2*pi*1.0;       % wavenumber (increase/decrease to change oscillation)
radius = 1.0;       % radius of sampling curve
cluster = 0.02;     % clustering strength: 0 (uniform) .. 0.1 (very clustered)
self_reg = 1e-1;    % diagonal regularization magnitude (physically from principal value)
add_noise = 0;      % small random symmetric perturbation (0..small)
make_harder = false; % if true: produce very close pairs to make it nastier

% build matrix
Z = build_EFIE_like(N, k, radius, cluster, self_reg, add_noise, make_harder);
ht = Z(:, 1);
%}
b = ones(N, 1);

tolerance = 1e-12;
maxIterCount = 1000;

%{
% Fourier coeffs of x^2 on [0, 2pi]
ks = 0:(N-1);              % nonnegative k for first column
ht = zeros(N,1);            % will hold t_0, t_1, ..., t_{N-1}

% numerical integration: vectorized loop with integral
for idx = 1:numel(ks)
    k = ks(idx);
    if k == 0
        ht(idx) = (1/(2*pi)) * integral(@(x) x.^4, 0, 2*pi);  % real scalar
    else
        integrand = @(x) (x.^4).*exp(-1i * k * x);
        ht(idx) = (1/(2*pi)) * integral(integrand, 0, 2*pi);
    end
end

% Right-hand side vector
b = ones(N, 1);

T = toeplitz(ht, ht);

%% mlin direct solver
tic;
x1 = T \ b;
toc;

%% Plain CG solver
tolerance = 1e-12;
maxIterCount = 100;

tic;
[x2, flagC, relresC, iterC] = cgs(T, b, tolerance, maxIterCount);
tcgs = toc;
fprintf('Plain cgs: flag=%d, relres=%g, iter=%d, time=%g s\n', flagC, relresC, iterC, tcgs);

% Check the solution distance.
fprintf('||x_direct - x_cgs|| = %g\n', norm(x1 - x2));

%% Strang's preconditioner
S = zeros(N, N);
s = zeros(N, 1);

half = floor(N / 2);

for n = 0:(N - 1)
    if n <= half
        s(n + 1) = ht(n + 1);
    else
        s(n + 1) = ht(N - n + 1);
    end
end

xS = zeros(N, 1);
tic;
%[xS, flagC, relresC, iterC] = cgs(T, b, tolerance, maxIterCount, S);
tcgs = toc;
fprintf('Strang-preconditioned cgs: flag=%d, relres=%g, iter=%d, time=%g s\n', flagC, relresC, iterC, tcgs);

% Check the solution distance.
fprintf('||x_direct - x_strang_cgs|| = %g\n', norm(x1 - xS));

%% T. Chan's preconditioner
% Build the first column c of the circulant c(A) using formula:
% c_j = ((n-j)*a_j + j*a_{j-n})/n  for 0 <= j < n
% with a_{j-n} == a_{n-j} for symmetric Toeplitz.
c = zeros(N,1);
for j = 0:(N-1)
    if j == 0
        c(1) = ht(1);   % j = 0 => c0 = a0
    else
        % a_j is a(j+1); a_{j-n} = a_{n-j} = a(N-j+1)
        c(j+1) = ((N - j) * ht(j+1) + j * ht(N - j + 1)) / N;
    end
end

xC = zeros(N, 1);
tic;
%[xC, flagC, relresC, iterC] = cgs(T, b, tolerance, maxIterCount, C);
tcgs = toc;
fprintf('Chan-preconditioned cgs: flag=%d, relres=%g, iter=%d, time=%g s\n', flagC, relresC, iterC, tcgs);

% Compare solutions (to direct solve)
fprintf('||x_direct - x_chan_cgs|| = %g\n', norm(x1 - xC));
%}

s = zeros(N, 1);
half = floor(N / 2);
for n = 0:(N - 1)
    if n <= half
        s(n + 1) = ht(n + 1);
    else
        s(n + 1) = ht(N - n + 1);
    end
end

c = zeros(N,1);
for j = 0:(N-1)
    if j == 0
        c(1) = ht(1);   % j = 0 => c0 = a0
    else
        % a_j is a(j+1); a_{j-n} = a_{n-j} = a(N-j+1)
        c(j+1) = ((N - j) * ht(j+1) + j * ht(N - j + 1)) / N;
    end
end
%{
%% Plain CG solver with FFT mat-vec multiplication
data = FFTMultPrep(ht);
applicator = @(x) FFTMult(data, x);

% Eigenvalues
eigsFun = @(x) FFTMult(data, x);
eigs_T = eigs(eigsFun, N)

tic;
[x3, flagC, relresC, iterC] = cgs(applicator, b, tolerance, maxIterCount);
tcgs = toc;
fprintf('FFT accelerated plain cgs: flag=%d, relres=%g, iter=%d, time=%g s\n', flagC, relresC, iterC, tcgs);

% Compare solutions (to direct solve)
fprintf('||x_3 - x_fft_cgs|| = %g\n', norm(x3 - x3));

%% Strang-preconditioned CG solver with FFT mat-vec multiplication
sf = fft(s);

function y = applyPreconditioner(sf, x)
    y = ifft(fft(x) ./ sf);
end

precondApplicator = @(x) applyPreconditioner(sf, x);

tic;
[xSf, flagC, relresC, iterC] = cgs(applicator, b, tolerance, maxIterCount, precondApplicator);
tcgs = toc;
fprintf('FFT accelerated Strang-preconditioned cgs: flag=%d, relres=%g, iter=%d, time=%g s\n', flagC, relresC, iterC, tcgs);

% Compare solutions (to direct solve)
fprintf('||x_3 - x_fft_strang_cgs|| = %g\n', norm(x3 - xSf));

% Eigenvalues
function y = eigsFunStrang_(data, sf, x)
    y = applyPreconditioner(sf, FFTMult(data, x));
end

eigsFunStrang = @(x) eigsFunStrang_(data, sf, x);

eigs_Strang = eigs(eigsFunStrang, N)

%% T. Chan-preconditioned CG solver with FFT mat-vec multiplication
cf = fft(c);

precondApplicator = @(x) applyPreconditioner(cf, x);

tic;
[xCf, flagC, relresC, iterC] = cgs(applicator, b, tolerance, maxIterCount, precondApplicator);
tcgs = toc;
fprintf('FFT accelerated Chan-preconditioned cgs: flag=%d, relres=%g, iter=%d, time=%g s\n', flagC, relresC, iterC, tcgs);

% Compare solutions (to direct solve)
fprintf('||x_3 - x_fft_chan_cgs|| = %g\n', norm(x3 - xCf));

% Eigenvalues
function y = eigsFunChan_(data, cf, x)
    y = applyPreconditioner(cf, FFTMult(data, x));
end

eigsFunChan = @(x) eigsFunChan_(data, cf, x);

eigs_Chan = eigs(eigsFunChan, N)
%}
%% Plain CG solver with FFT mat-vec multiplication
data = FFTMultPrep(ht);
applicator = @(x) FFTMult(data, x);

function y = applyPreconditioner(sf, x)
    y = ifft(fft(x) ./ sf);
end

%% New preconditioner: T + S is a good preconditioner
% ht is the first column of the symmetric Toeplitz matrix
st = [0; flipud(ht(2:end))];
N1 = ht + st;
N1f = fft(N1);

precondApplicator = @(x) applyPreconditioner(N1f, x);

tic;
[xNew1, flagC, relresC, iterC] = cgs(applicator, b, tolerance, maxIterCount, precondApplicator);
tcgs = toc;
fprintf('FFT accelerated new1-preconditioned cgs: flag=%d, relres=%g, iter=%d, time=%g s\n', flagC, relresC, iterC, tcgs);

% Eigenvalues
function y = eigsFunNew1_(data, N1f, x)
    y = applyPreconditioner(N1f, FFTMult(data, x));
end

eigsFunNew1 = @(x) eigsFunNew1_(data, N1f, x);

eigs_New1 = eigs(eigsFunNew1, N)

%% New preconditioner: T - S is a good preconditioner
% ht is the first column of the symmetric Toeplitz matrix
C = zeros(N, N);
for n = 0:(N - 1)
    C(:, n + 1) = circshift(c, n);
end
C

r1 = C * b;

k = (0:(N - 1)).';
d = exp(1i * pi * k / N)

C = diag(d) * C * diag(conj(d))

data.hsf = fft(hs);
data.d = d;

r2 = FFTSkewMult(data, b)

function y = FFTSkewMult(data, x)
    d = data.d;
    dx = d .* x;
    y = conj(d) .* ifft(data.hsf .* fft(dx));
end


precondApplicator = @(x) applyPreconditioner(N1f, x);

tic;
[xNew1, flagC, relresC, iterC] = cgs(applicator, b, tolerance, maxIterCount, precondApplicator);
tcgs = toc;
fprintf('FFT accelerated new1-preconditioned cgs: flag=%d, relres=%g, iter=%d, time=%g s\n', flagC, relresC, iterC, tcgs);

% Eigenvalues
function y = eigsFunNew2_(data, N1f, x)
    y = applyPreconditioner(N1f, FFTMult(data, x));
end

eigsFunNew1 = @(x) eigsFunNew2_(data, N1f, x);

eigs_New1 = eigs(eigsFunNew1, N)
%{
%% New preconditioner: C_1^(-1) is a good preconditioner
function y = applyPreconditioner2(data, x)
    N = data.N;
    y = ifft(fft(x, 2*N) ./ data.hcf);
    y = y(1:N);
end

precondApplicator2 = @(x) applyPreconditioner2(data, x);

tic;
[xNew2, flagC, relresC, iterC] = cgs(applicator, b, tolerance, maxIterCount, precondApplicator2);
tcgs = toc;
fprintf('FFT accelerated new2-preconditioned cgs: flag=%d, relres=%g, iter=%d, time=%g s\n', flagC, relresC, iterC, tcgs);

% Compare solutions (to direct solve)
fprintf('||x_3 - x_fft_new2_cgs|| = %g\n', norm(x3 - xNew2));

% Eigenvalues
function y = eigsFunNew2_(data, x)
    y = applyPreconditioner2(data, FFTMult(data, x));
end

eigsFunNew2 = @(x) eigsFunNew2_(data, x);

eigs_New2 = eigs(eigsFunNew2, N)

%% New preconditioner: full 
function y = applyPreconditioner3(data, x)
    N = data.N;
    % 1) u = C^{-1} x  (via FFT on length 2N)
    u_big = ifft( fft(x,2*N) ./ data.hcf );
    u = u_big(1:N);

    % 2) w = T * u  (Toeplitz multiply)
    w = FFTMult(data, u);

    w = 2 * x - w;

    % 3) z = C^{-1} w
    y_big = ifft( fft(w,2*N) ./ data.hcf );
    y = y_big(1:N);
end

precondApplicator3 = @(x) applyPreconditioner3(data, x);

tic;
[xNew3, flagC, relresC, iterC] = cgs(applicator, b, tolerance, maxIterCount, precondApplicator3);
tcgs = toc;
fprintf('FFT accelerated new3-preconditioned cgs: flag=%d, relres=%g, iter=%d, time=%g s\n', flagC, relresC, iterC, tcgs);

% Compare solutions (to direct solve)
fprintf('||x_3 - x_fft_new3_cgs|| = %g\n', norm(x3 - xNew3));

% Eigenvalues
function y = eigsFunNew3_(data, x)
    y = applyPreconditioner3(data, FFTMult(data, x));
end

eigsFunNew3 = @(x) eigsFunNew3_(data, x);

eigs_New3 = eigs(eigsFunNew3, N)

%% Nearest circulant preconditioner 
% ht : first column of symmetric Toeplitz (rho_0, rho_1, ..., rho_{N-1})
% Example: ht = [2; 0.8; 0.3; 0.1];  % replace with your data

ht = ht(:);            % ensure column
N = numel(ht);

% build c (first column of optimal circulant C_Sigma)
c = zeros(N,1);
for m = 0:(N-1)
    if m == 0
        c(1) = ht(1);                    % c0 = rho0
    else
        % indices: rho_m = ht(m+1), rho_{N-m} = ht(N-m+1)
        c(m+1) = ht(m+1) + (m / N) * ( ht(N - m + 1) - ht(m+1) );
    end
end

cf = fft(c);

precondApplicator = @(x) applyPreconditioner(cf, x);

tic;
[xCf, flagC, relresC, iterC] = cgs(applicator, b, tolerance, maxIterCount, precondApplicator);
tcgs = toc;
fprintf('FFT accelerated Nearest-circulant_preconditioned cgs: flag=%d, relres=%g, iter=%d, time=%g s\n', flagC, relresC, iterC, tcgs);

% Compare solutions (to direct solve)
fprintf('||x_3 - x_fft_nearest_circulant_cgs|| = %g\n', norm(x3 - xCf));
%}
%% Helpers
function data = FFTMultPrep(ht)
    data.N = numel(ht);
    data.hcf = fft([ht; 0; flipud(ht(2:end))]);
end

% 1D fft accelerated muliplication
function y = FFTMult(data, x)
    N = data.N;
    y = ifft(data.hcf .* fft(x, 2*N));
    y = y(1:N);
end

function Z = build_EFIE_like(N, k, radius, cluster_strength, self_reg, add_noise, make_harder)
% build N sample points on circle, optionally cluster them
theta = linspace(0, 2*pi, N+1).';
theta(end) = [];
% introduce small random jitter; clustering makes some nearby
if cluster_strength > 0
    jitter = cluster_strength * (0.5 - rand(N,1)) .* (1 + 10*rand(N,1).*(0.1*make_harder));
else
    jitter = zeros(N,1);
end
theta = theta + jitter;

% optionally create very close pairs to make conditioning worse
if make_harder
    % take a random small set and cluster them tightly
    m = max(2, round(0.02*N)); % number of tight cluster points
    idx = randperm(N, m);
    base = theta(idx(1));
    theta(idx) = base + (1e-4)*(rand(m,1)-0.5); % extremely close
end

x = radius * cos(theta);
y = radius * sin(theta);
pts = [x, y];

% build distance matrix
DX = x(:) - x(:).';
DY = y(:) - y(:).';
R = sqrt(DX.^2 + DY.^2);

% kernel (2D/3D free-space choice; simple scalar Green's function)
% note: avoid dividing by zero on diagonal
eps_replace = 1e-6;
Rdiag = max(R, eps_replace);

G = exp(-1i * k * Rdiag) ./ Rdiag;

% regularize diagonal (approx self-impedance / principal value)
% set diagonal to a physically plausible complex self-term
for ii=1:N
    G(ii,ii) = (1 - 1i) * self_reg;   % tune self_reg to change conditioning
end

% ensure complex-symmetric (EFIE-like)
Z = 0.5*(G + G.');

% optionally add small symmetric noise (modeling discretization error)
if add_noise > 0
    E = add_noise * (randn(N) + 1i*randn(N));
    E = 0.5*(E + E.'); % symmetric
    Z = Z + E;
end
end
