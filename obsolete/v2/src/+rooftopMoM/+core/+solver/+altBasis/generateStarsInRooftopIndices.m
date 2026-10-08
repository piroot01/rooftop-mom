function starsConverted = generateStarsInRooftopIndices(countX, countY)
% generateStarsInRooftopIndices: Returns 4xN matrix describing which rooftops
%   makes erach star basis function.
%
%   Inputs:
%       countX: Positive integer specifying the number of patches in the x direction.
%       countY: Positive integer specifying the number of patches in the y direction.
%
%   Output:
%       stars: 4xN matrix describing the star basis functions.

    basisSizeX = (countX - 1) * countY;
    basisSizeY = countX * (countY - 1);

    % Create a matrix of the sequential indices for the vertical edges.
    xIndices = reshape(1:basisSizeX, countX - 1, countY);
    rightMat = [xIndices; zeros(1, countY)];
    right = rightMat(:)';

    leftMat = [zeros(1, countY); xIndices];
    left = leftMat(:)';

    yIndices = basisSizeX + (1:basisSizeY);

    top = [yIndices, zeros(1, countX)];

    bottom = [zeros(1, countX), yIndices];

    starsConverted = [right; top; left; bottom];
end
