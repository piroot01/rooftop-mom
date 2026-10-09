function [blockXX, blockXY, blockYX, blockYY] = generateSpatialBlocks(generators, faceSizeX, faceSizeY)
% generateSpatialBlocks: Returns the two-sided generators of the four 2x2 blocks of the Green's matrix on the face grid.
%   The scalar potential enters through differences of the electric generator e; at the offset (p, q) between the first
%   faces of two rooftops
%       XX(p, q) = magneticXX(p, q) - (2 * e(p, q) - e(p - 1, q) - e(p + 1, q)) / faceSizeX^2,
%       XY(p, q) = -(e(p, q) - e(p, q - 1) - e(p + 1, q) + e(p + 1, q - 1)) / (faceSizeX * faceSizeY),
%   YY follows as XX along y, and YX(p, q) = XY(-p, -q) as the matrix is symmetric. The offsets that no pair of rooftops
%   of a block has are zero.
%
%   Inputs:
%       generators: rooftopMoM.assembler.greensMatrix.Generators.
%       faceSizeX: Face size in the x direction.
%       faceSizeY: Face size in the y direction.
%
%   Outputs:
%       blockXX, blockXY, blockYX, blockYY: (2 * countX - 1) x (2 * countY - 1) two-sided generators, empty for a
%           missing direction.

    electric = generators.electric();
    [countX, countY] = size(electric);
    hasX = generators.content() ~= rooftopMoM.basisFunctions.RooftopContent.YONLY;
    hasY = generators.content() ~= rooftopMoM.basisFunctions.RooftopContent.XONLY;

    [blockXX, blockXY, blockYX, blockYY] = deal([]);

    % A border of zeros lets e be shifted by one face; shifted(a, b) holds e(p + a, q + b).
    bordered = zeros(2 * countX + 1, 2 * countY + 1, 'like', electric);
    bordered(2:(end - 1), 2:(end - 1)) = rooftopMoM.assembler.preconditioner.utils.extendEvenly(electric);
    shifted = @(a, b) bordered((2:(2 * countX)) + a, (2:(2 * countY)) + b);

    % The offsets p along the rows and q along the columns.
    p = (-(countX - 1):(countX - 1)).';
    q = -(countY - 1):(countY - 1);

    if hasX
        blockXX = extendMagnetic_(generators.magneticXX(), countX, countY) ...
            - (2 * shifted(0, 0) - shifted(-1, 0) - shifted(1, 0)) / faceSizeX^2;
        blockXX(abs(p) > countX - 2, :) = 0;
    end

    if hasY
        blockYY = extendMagnetic_(generators.magneticYY(), countX, countY) ...
            - (2 * shifted(0, 0) - shifted(0, -1) - shifted(0, 1)) / faceSizeY^2;
        blockYY(:, abs(q) > countY - 2) = 0;
    end

    if hasX && hasY
        blockXY = -(shifted(0, 0) - shifted(0, -1) - shifted(1, 0) + shifted(1, -1)) / (faceSizeX * faceSizeY);

        % An x-directed rooftop is at most countX - 2 faces right of a y-directed one, which is at most countY - 2
        % faces above it.
        blockXY(p > countX - 2, :) = 0;
        blockXY(:, q < -(countY - 2)) = 0;
        blockYX = rot90(blockXY, 2);
    end
end

function twoSided = extendMagnetic_(magnetic, countX, countY)
    twoSided = rooftopMoM.assembler.preconditioner.utils.extendEvenly( ...
        rooftopMoM.assembler.utils.padWithZeros(magnetic, countX, countY));
end
