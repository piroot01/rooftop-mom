function circulant = generateRChan(generator)
% generateRChan: Returns the R. Chan circulant approximation of a symmetric Toeplitz generator: the generator plus its
%   point reflection (offsets counted from the end) in the part with both offsets nonzero.
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

    reflection = zeros(size(generator), 'like', generator);
    reflection(2:end, 2:end) = rot90(generator(2:end, 2:end), 2);

    circulant = generator + reflection;
end
