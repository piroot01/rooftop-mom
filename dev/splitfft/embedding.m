clc; clearvars;

N = 3;

x = randi(10, [N, 1]);

ht1 = randi(10, [N, 1]);
ht2 = randi(10, [N, 1]);
ht1(1) = ht2(1);

ht1 = ht2;

ht = [ht2(end:(-1):2); ht1];

T = toeplitz(ht1, ht2)

T * x

U = diag(ones(2*N - 2,1 ), 1);

P = diag(ones(2*N - 2, 1), -1);
P(1, end) = 1;

ht;

P * ht;


Cht1 = circulant(ht1, 1)
C = circulant(ht, 1)

B = ifft(fft(circshift(ht, 0)) .* fft([x; zeros(N - 1, 1)]))

b = B(N:end)


x;

x .* [0, 0, 1];

P

R = [zeros(N, N - 1), eye(N)];

L = [eye(N), zeros(N, N - 1)];

h0 = circshift(ht, N);

L * U^1 * h0

L * U^1

L * P^1

%R * P^0

%R * P^1

%R * P^2

%%
clc;

N = 4;

h0 = 1:N;

T = toeplitz(h0)

h = [h0, h0(end:(-1):2)];
h = h(:);

U = diag(ones(2*N - 2,1 ), 1);

P = diag(ones(2*N - 2, 1), -1);
P(1, end) = 1;

R = [zeros(N, N - 1), eye(N)];

L = [eye(N), zeros(N, N - 1)];
%%
clc;

N = 4;
k = 100;
h0 = 1:N;
h = [h0, zeros(1, k), h0(end:(-1):2)];
h = h(:);

T = toeplitz(h0);

x = randi(10, [N, 1]);

T * x

K = numel(h);

P = diag(ones(K - 1, 1), -1);
P(1, end) = 1;

L = [eye(N), zeros(N, K - N)];

ifft(fft(h) .* fft([x; zeros([K - N, 1])]))

