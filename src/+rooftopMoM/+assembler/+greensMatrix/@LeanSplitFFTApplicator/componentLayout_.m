function [rowCounts, columnCounts, offsets] = componentLayout_(obj)
% componentLayout_: Returns the grid layout of the x- and y-directed rooftop components.
%
%   Outputs:
%       rowCounts: 1x2 vector of the grid rows of the x- and y-directed components.
%       columnCounts: 1x2 vector of the grid columns of the x- and y-directed components.
%       offsets: 1x2 vector of the offsets of the components in the rooftop vectors.

    [rowCountX, columnCountX] = obj.rooftops_.gridSize('x');
    [rowCountY, columnCountY] = obj.rooftops_.gridSize('y');

    rowCounts = [rowCountX, rowCountY];
    columnCounts = [columnCountX, columnCountY];
    offsets = [0, rowCountX * columnCountX];
end
