clear; clearvars;

%rng(0);

% Construction of the 4d Toeplitz tensor kernel.
% Dimensions of the kernel.
M = 100;
N = 2;

repeatCount = 1000;

ht = randi(100, [M, N]) + 1i * randi(100, [M, N]);

% Embedding to the circulant tensor kernel.
rIdx = [1:M, (M - 1):(-1):2];
cIdx = [1:N, (N - 1):(-1):2];

hc0 = ht(rIdx, cIdx);

% Convert the circulant tensor to spectral domain.
hc_fft = fft2(hc0);

% Construct the vector that will be multiplied.
x = randi(100, [M*N, 1]);

%% 1. method: Classical vector embedding and reshaping
% Reshape the vector
t0 = tic;

for n = 1:repeatCount

x_reshaped = reshape(x, M, N);

% Padd and apply fft.
x_reshaped_fft = fft2(x_reshaped, 2*M - 2, 2*N - 2);

b_raw = ifft2(hc_fft .* x_reshaped_fft);
b_cutted = b_raw(1:M, 1:N);

b0 = b_cutted(:);

end

time0 = toc(t0);

%% 2. method: Split fft
% Preparation of the circulant kernel
% Embed the Toeplitz kernel into the circulant one with a size 2N x 2M
% ht is M-by-N
prep = splitFFTMultDataPrep(ht);

t1 = tic;

for n = 1:repeatCount
    b1 = splitFFTMult(prep, x);
end

time1 = toc(t1);

max(abs(b0 - b1))

fprintf('time0 = %g s\ntime1 = %g s\n', time0 / repeatCount, time1 / repeatCount);
fprintf('%gx speed up over classical implementation\n', round(time0 / time1, 2));
