function result = originalFFTMult(data, input)
    M = data.M; N = data.N;

    inputReshaped = reshape(input, M, N);

    inputReshapedFFT = fft2(inputReshaped, 2 * M - 2, 2 * N - 2);

    result = ifft2(data.circulantKernelFFT .* inputReshapedFFT);

    result = result(1:M, 1:N);

    result = result(:);
end