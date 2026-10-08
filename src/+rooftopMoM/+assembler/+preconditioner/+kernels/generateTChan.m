function circulant = generateTChan(generator)
% generateTChan: Returns the optimal circulant approximation of T. Chan of a symmetric Toeplitz generator, the
%   circulant closest in the Frobenius norm. Each offset averages the generator with its wrapped counterpart,
%   weighted by how often each occurs in the matrix.
%
%   Input:
%       generator: 2D Toeplitz generator (may be empty).
%
%   Output:
%       circulant: Generator of the circulant, of the size of the input.

    if isempty(generator)
        circulant = [];
        return;
    end

    [countX, countY] = size(generator);
    p = (0:(countX - 1)).';
    q = 0:(countY - 1);
    wrappedP = mod(countX - p, countX) + 1;
    wrappedQ = mod(countY - q, countY) + 1;

    circulant = ((countX - p) .* (countY - q) .* generator(p + 1, q + 1) ...
        + p .* (countY - q) .* generator(wrappedP, q + 1) ...
        + (countX - p) .* q .* generator(p + 1, wrappedQ) ...
        + p .* q .* generator(wrappedP, wrappedQ)) / (countX * countY);
end
