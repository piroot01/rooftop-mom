clc; clearvars;

N = 16;

ht = zeros(N, 1);
ht(1) = 2;
ht(2) = -1;

T = toeplitz(ht, ht);
Teig = eig(T);


hc = 
function fftMult