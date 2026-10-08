function [blockXX, blockXY, blockYX, blockYY] = generateSpatialBlocks(generators, faceSizeX, faceSizeY)
% generateSpatialBlocks: Returns the Toeplitz generators of the four 2x2 blocks of the Green's matrix on the face grid.
%   The scalar potential enters through differences of the electric generator; its even extension closes the
%   stencils at the first offset, a reflection at the last offset closes them there.
%
%   Inputs:
%       generators: rooftopMoM.assembler.greensMatrix.Generators.
%       faceSizeX: Face size in the x direction.
%       faceSizeY: Face size in the y direction.
%
%   Outputs:
%       blockXX, blockXY, blockYX, blockYY: countX x countY generators, empty for a missing direction.

    electric = generators.electric();
    [countX, countY] = size(electric);
    hasX = generators.content() ~= rooftopMoM.basisFunctions.RooftopContent.YONLY;
    hasY = generators.content() ~= rooftopMoM.basisFunctions.RooftopContent.XONLY;

    blockXX = [];
    blockXY = [];
    blockYX = [];
    blockYY = [];

    if hasX
        difference = zeros(countX, countY, 'like', electric);
        difference(2:(end - 1), :) = 2 * electric(2:(end - 1), :) - electric(1:(end - 2), :) - electric(3:end, :);
        difference(1, :) = 2 * (electric(1, :) - electric(2, :));
        difference(end, :) = 2 * (electric(end, :) - electric(end - 1, :));

        blockXX = rooftopMoM.assembler.utils.padWithZeros(generators.magneticXX(), countX, countY) ...
            - (1 / faceSizeX)^2 .* difference;
    end

    if hasY
        difference = zeros(countX, countY, 'like', electric);
        difference(:, 2:(end - 1)) = 2 * electric(:, 2:(end - 1)) - electric(:, 1:(end - 2)) - electric(:, 3:end);
        difference(:, 1) = 2 * (electric(:, 1) - electric(:, 2));
        difference(:, end) = 2 * (electric(:, end) - electric(:, end - 1));

        blockYY = rooftopMoM.assembler.utils.padWithZeros(generators.magneticYY(), countX, countY) ...
            - (1 / faceSizeY)^2 .* difference;
    end

    if hasX && hasY
        % XY: the difference E(i, j - 1) - E(i, j) in y, then its difference in x.
        differenceY = zeros(countX, countY, 'like', electric);
        differenceY(:, 2:end) = electric(:, 1:(end - 1)) - electric(:, 2:end);
        differenceY(:, 1) = -differenceY(:, 2);

        mixed = zeros(countX, countY, 'like', electric);
        mixed(1:(end - 1), :) = differenceY(1:(end - 1), :) - differenceY(2:end, :);
        mixed(end, :) = -mixed(end - 1, :);
        blockXY = (1 / (faceSizeX * faceSizeY)) .* mixed;

        % YX: the difference E(i - 1, j) - E(i, j) in x, then its difference in y.
        differenceX = zeros(countX, countY, 'like', electric);
        differenceX(2:end, :) = electric(1:(end - 1), :) - electric(2:end, :);
        differenceX(1, :) = -differenceX(2, :);

        mixed = zeros(countX, countY, 'like', electric);
        mixed(:, 1:(end - 1)) = differenceX(:, 1:(end - 1)) - differenceX(:, 2:end);
        mixed(:, end) = -mixed(:, end - 1);
        blockYX = (1 / (faceSizeY * faceSizeX)) .* mixed;
    end
end
