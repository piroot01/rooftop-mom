function data = originalFFTMultDataPrep(toeplitzKernel)
    [M, N] = size(toeplitzKernel);
    data.M = M; data.N = N;

    rowIdx = [1:M, (M - 1):(-1):2];
    columnIdx = [1:N, (N - 1):(-1):2];

    data.circulantKernelFFT = fft2(toeplitzKernel(rowIdx, columnIdx));
end