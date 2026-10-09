function eigenvalues = generateOptimalOfProduct(partsA, partsB)
% generateOptimalOfProduct: Returns the eigenvalues of c(A * B'), the optimal circulant of T. Chan of the product of
%   two two-level Toeplitz matrices, from the eigenvalues of their parts (generatePartSpectra), without the product.
%   Along each direction c keeps a circulant factor, c(C * M) = C * c(M), and projects a skew-circulant one by the
%   circular convolution of its eigenvalues with
%       K(d) = 1 / (count^2 * sin(pi * (2 * d + 1) / (2 * count))^2),  d = 0, ..., count - 1.
%   So every pair of parts adds the eigenvalues of the A part times the conjugate eigenvalues of the B part. Along a
%   direction where one part is skew-circulant, its eigenvalues are convolved before the product; where both are,
%   the product is convolved.
%
%   Inputs:
%       partsA: 2x2 cell of the eigenvalues of the parts of A.
%       partsB: 2x2 cell of the eigenvalues of the parts of B.
%
%   Output:
%       eigenvalues: countX x countY eigenvalues of c(A * B').

    [countX, countY] = size(partsA{1, 1});
    kernels = {kernel_(countX), kernel_(countY)};

    % The parts along x and y: 1 for circulant, 2 for skew-circulant.
    types = [1, 1; 2, 1; 1, 2; 2, 2];
    eigenvalues = zeros(countX, countY, 'like', 1i);

    for m = 1:4
        for n = 1:4
            typeA = types(m, :);
            typeB = types(n, :);
            valuesA = partsA{typeA(1), typeA(2)};
            valuesB = partsB{typeB(1), typeB(2)};

            for dimension = find(typeA == 2 & typeB == 1)
                valuesA = convolve_(valuesA, kernels{dimension}, dimension);
            end

            for dimension = find(typeA == 1 & typeB == 2)
                valuesB = convolve_(valuesB, kernels{dimension}, dimension);
            end

            product = valuesA .* conj(valuesB);

            for dimension = find(typeA == 2 & typeB == 2)
                product = convolve_(product, kernels{dimension}, dimension);
            end

            eigenvalues = eigenvalues + product;
        end
    end
end

function spectrum = kernel_(count)
    distances = (0:(count - 1)).';
    spectrum = fft(1 ./ (count^2 * sin(pi * (2 * distances + 1) / (2 * count)).^2));
end

function values = convolve_(values, kernel, dimension)
    if dimension == 1
        values = ifft(kernel .* fft(values, [], 1), [], 1);
    else
        values = ifft(kernel.' .* fft(values, [], 2), [], 2);
    end
end
