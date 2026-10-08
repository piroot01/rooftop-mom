function shifts = generateRectGridShifts(countX, countY, sizeX, sizeY)
% generateRectGridShifts: Generates shift vectors for a rectangular grid.
%
%   Inputs:
%       countX: Number of cells in the x-direction.
%       countY: Number of cells in the y-direction.
%       sizeX:  Width of each cell in the x-direction.
%       sizeY:  Height of each cell in the y-direction.
%
%   Outputs:
%       shifts: A 3xN matrix where each column is the [x; y; 0] coordinate of the bottom-left corner of a cell.

    [xn, yn] = ndgrid(double(0:(countX - 1)) * sizeX, double(0:(countY - 1)) * sizeY);

    N = numel(xn);
    shifts = [reshape(xn, 1, N); reshape(yn, 1, N); zeros(1, N)];
end
