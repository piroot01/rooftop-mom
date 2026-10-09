function [rowCounts, columnCounts, offsets] = componentLayout_(rooftops)
% componentLayout_: Returns the grid layout of the x- and y-directed rooftop components.
%
%   Input:
%       rooftops: rooftopMoM.basisFunctions.UniformRooftops.
%
%   Outputs:
%       rowCounts: 1x2 vector of the grid rows of the x- and y-directed components.
%       columnCounts: 1x2 vector of the grid columns of the x- and y-directed components.
%       offsets: 1x2 vector of the offsets of the components in the rooftop vectors.

    [rowCountX, columnCountX] = rooftops.gridSize('x');
    [rowCountY, columnCountY] = rooftops.gridSize('y');

    rowCounts = [rowCountX, rowCountY];
    columnCounts = [columnCountX, columnCountY];
    offsets = [0, rowCountX * columnCountX];
end
