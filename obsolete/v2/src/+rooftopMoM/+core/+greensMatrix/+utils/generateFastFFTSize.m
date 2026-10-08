function fftSize = generateFastFFTSize(count)
% generateFastFFTSize: Returns the smallest FFT size not below count whose prime factors are at most 7.
%
%   Input:
%       count: Positive integer specifying the minimal required size.
%
%   Output:
%       fftSize: Positive integer, count rounded up to a 7-smooth number (FFTW is fastest for such sizes).

    fftSize = count;

    while max(factor(fftSize)) > 7
        fftSize = fftSize + 1;
    end
end
