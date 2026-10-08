function circulant = generateRChan(generator)
% generateRChan: Returns the R. Chan circulant approximation of a symmetric Toeplitz generator: in every direction the
%   offset k > 0 gets the generator at k plus the generator at count - k, so the circulant stays symmetric.
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

    circulant = addReflection_(addReflection_(generator, 1), 2);
end

function reflected = addReflection_(generator, dimension)
    count = size(generator, dimension);
    indices = repmat({':'}, 1, 2);
    indices{dimension} = [1, count:(-1):2];
    reflection = generator(indices{:});

    % The zero offset has no counterpart.
    indices{dimension} = 1;
    reflection(indices{:}) = 0;

    reflected = generator + reflection;
end
