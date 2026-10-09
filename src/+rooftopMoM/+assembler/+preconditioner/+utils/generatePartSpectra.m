function spectra = generatePartSpectra(generator)
% generatePartSpectra: Returns the eigenvalues of the four parts of a two-level Toeplitz matrix. Along each direction a
%   Toeplitz matrix is the sum of a circulant and a skew-circulant matrix, T = C + S, with the generators
%       c(k) = (t(k) + t(k - count)) / 2,  s(k) = (t(k) - t(k - count)) / 2,  k = 0, ..., count - 1,
%   where the missing t(-count) is zero. A circulant is diagonal in the Fourier basis, a skew-circulant in the Fourier
%   basis modulated by exp(1i * pi * k / count).
%
%   Input:
%       generator: (2 * countX - 1) x (2 * countY - 1) two-sided generator.
%
%   Output:
%       spectra: 2x2 cell of countX x countY eigenvalues. The first index is the part along x, the second along y:
%           1 for the circulant and 2 for the skew-circulant part.

    countX = (size(generator, 1) + 1) / 2;
    countY = (size(generator, 2) + 1) / 2;
    modulationX = exp(1i * pi * (0:(countX - 1)).' / countX);
    modulationY = exp(1i * pi * (0:(countY - 1)) / countY);

    partsX = splitRows_(generator);
    spectra = cell(2);

    for a = 1:2
        % The transpose splits along y, the columns.
        partsY = splitRows_(partsX{a}.');

        for b = 1:2
            part = partsY{b}.';

            if a == 2
                part = modulationX .* part;
            end

            if b == 2
                part = part .* modulationY;
            end

            spectra{a, b} = fft2(part);
        end
    end
end

function parts = splitRows_(generator)
    % The row k + 1 of a part holds its generator at the offset k.
    count = (size(generator, 1) + 1) / 2;
    positive = generator(count:end, :);
    negative = [zeros(1, size(generator, 2), 'like', generator); generator(1:(count - 1), :)];
    parts = {(positive + negative) / 2, (positive - negative) / 2};
end
