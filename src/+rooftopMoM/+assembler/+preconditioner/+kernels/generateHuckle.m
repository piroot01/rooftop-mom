function circulant = generateHuckle(generator, bandwidthX, bandwidthY)
% generateHuckle: Returns the Huckle circulant approximation of a symmetric Toeplitz generator: the offsets below the
%   bandwidth, damped linearly by (1 - offset / bandwidth) in every direction and placed at the offset and at its
%   wrapped counterpart.
%
%   Inputs:
%       generator: 2D Toeplitz generator (may be empty).
%       bandwidthX: Positive integer specifying the bandwidth in the x direction.
%       bandwidthY: Positive integer specifying the bandwidth in the y direction.
%
%   Output:
%       circulant: Generator of the circulant, of the size of the input.

    if isempty(generator)
        circulant = [];
        return;
    end

    [countX, countY] = size(generator);
    i = (0:min(bandwidthX - 1, countX - 1)).';
    j = 0:min(bandwidthY - 1, countY - 1);
    damped = (1 - i / bandwidthX) .* (1 - j / bandwidthY) .* generator(i + 1, j + 1);

    rowsX = {i + 1, mod(countX - i, countX) + 1};
    columnsY = {j + 1, mod(countY - j, countY) + 1};
    circulant = zeros(countX, countY, 'like', generator);

    for a = 1:2
        for b = 1:2
            % Offset 0 has no wrapped counterpart, so it is placed once.
            keepI = a == 1 | i > 0;
            keepJ = b == 1 | j > 0;
            rows = rowsX{a}(keepI);
            columns = columnsY{b}(keepJ);
            circulant(rows, columns) = circulant(rows, columns) + damped(keepI, keepJ);
        end
    end
end
