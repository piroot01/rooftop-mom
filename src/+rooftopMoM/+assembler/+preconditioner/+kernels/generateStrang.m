function circulant = generateStrang(generator)
% generateStrang: Returns the Strang circulant approximation of a symmetric Toeplitz generator: the generator up to
%   half the size in every direction and its mirror image beyond.
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
    circulant = generator(mirror_(countX), mirror_(countY));
end

function indices = mirror_(count)
    offsets = 0:(count - 1);
    indices = min(offsets, count - offsets) + 1;
end
