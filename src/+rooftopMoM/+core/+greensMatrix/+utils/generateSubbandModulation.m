function modulation = generateSubbandModulation(fftSize, dim)
% generateSubbandModulation: Returns the modulation exp(-1i * pi * d / fftSize), d = 0, ..., fftSize - 1, which
%   maps the odd (half-bin shifted) frequencies of a circulant with the period 2 * fftSize onto integer ones.
%
%   Inputs:
%       fftSize: Positive integer specifying the subband size.
%       dim: Direction of the modulation, 1 for a column vector and 2 for a row vector.
%
%   Output:
%       modulation: fftSize x 1 (dim = 1) or 1 x fftSize (dim = 2) modulation vector.

    shape = [1, 1];
    shape(dim) = fftSize;

    modulation = reshape(exp(-1i * pi * (0:(fftSize - 1)) / fftSize), shape);
end
