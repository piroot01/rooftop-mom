function fftSize = generateFastFFTSize(count)
% generateFastFFTSize: Returns the smallest FFT size not below count whose prime factors are at most 7. FFTW is
%   fastest for such sizes.
%
%   Input:
%       count: Positive integer specifying the minimal size.
%
%   Output:
%       fftSize: Positive integer.

    fftSize = count;

    while max(factor(fftSize)) > 7
        fftSize = fftSize + 1;
    end
end
