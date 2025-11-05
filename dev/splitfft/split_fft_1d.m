clc; clearvars;

% Matrix size N
N = 100;

repeat = 1000;

% Toeplitz matrix 1st column - generator
ht = randi(100, [N, 1]);

x = randi(100, [N, 1]);

% exact multiplication
T = toeplitz(ht);

T * x

% Experimental split FFT
hcM = [ht; 0; ht(end:(-1):2)];

hcM_fft = fft(hcM, 2*N, 1);
hcM_fft_even = hcM_fft(1:2:end);
hcM_fft_odd = hcM_fft(2:2:end);

x_padded = [x; zeros(N, 1)];

x_padded_fft = fft(x_padded);
x_padded_fft(1:2:end);
x_padded_fft(2:2:end);

k = (0:(N - 1)).';
P = exp(-1i * pi * k ./ N);

P .* conj(P)

% this is correct
x_fft_even = fft(x, N, 1);
x_fft_odd = fft(P .* x, N, 1);

hcM_fft_x_even = hcM_fft_even .* x_fft_even;
hcM_fft_x_odd = hcM_fft_odd .* x_fft_odd;

b2_even = ifft(hcM_fft_x_even, N, 1);
b2_odd = conj(P) .* ifft(hcM_fft_x_odd, N, 1);

b2 = 0.5 * (b2_even + b2_odd)

