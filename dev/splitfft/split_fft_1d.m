clc; clearvars;

% Matrix size N
N = 10000;

repeat = 1000;

% Toeplitz matrix 1st column - generator
ht = randi(100, [N, 1]);

x = randi(100, [N, 1]);

%{
% Matlab's build-in mat-vec multiplication.
% Construct whole Toeplitz matrix.
T = toeplitz(ht);

t0 = tic;

for n = 1:repeat

b0 = T * x;

end

matlabMultTime = toc(t0);

fprintf("matlabMultTime = %g s\n", matlabMultTime / repeat);
%}

% Classical FFT multiplication method.
% Construct the circulant kernel.
hc = [ht; ht((end - 1):(-1):2)];
hc_fft = fft(hc, [], 1);



t1 = tic;

for n = 1:repeat

% Pad the x with zeros.
x_padded = [x; zeros(N - 2, 1)];

x_fft = fft(x_padded, [], 1);

b1_padded = ifft(hc_fft .* x_fft, [], 1);

b1 = real(b1_padded(1:N));

end

fftMultTime = toc(t1);
fprintf("fftMultTime = %g s\n", fftMultTime / repeat);

% Compare the vector norms.
%b0b1absNormDiff = abs(norm(b0) - norm(b1))

% Experimental split FFT
hcM = [ht; 0; ht(end:(-1):2)];

hcM_fft = fft(hcM, 2*N, 1);
hcM_fft_even = hcM_fft(1:2:end);
hcM_fft_odd = hcM_fft(2:2:end);

k = (0:(N - 1)).';
P = exp(1i * pi * k ./ N);

t2 = tic;

for n = 1:repeat

x_fft_even = fft(x, N, 1);
x_fft_odd = fft(P .* x, N, 1);
% TODO: Proof why this is happening.
x_fft_odd = circshift(x_fft_odd, -1);

hcM_fft_x_even = hcM_fft_even .* x_fft_even;
hcM_fft_x_odd = hcM_fft_odd .* x_fft_odd;

b2_even = ifft(hcM_fft_x_even, N, 1);
b2_odd = P .* ifft(hcM_fft_x_odd, N, 1);

b2 = 0.5 * (b2_even + b2_odd);

end

fftExpMultTime = toc(t2);

fprintf("fftExpMultTime = %g s\n", fftExpMultTime / repeat);
fprintf("%g speed increase", round(fftMultTime / fftExpMultTime, 2));

% Compare the vector norms.
b1b2absNormDiff = abs(norm(b1) - norm(b2))
