function shifts = generateRectGridShiftsSlice(countX, sizeX, sizeY, slice)
% generateRectGridShifts: Generates shift vectors for a rectangular grid.
%
%   Inputs:
%       countX: Number of cells in the x-direction.
%       sizeX:  Width of each cell in the x-direction.
%       sizeY:  Height of each cell in the y-direction.
%       slice: Indices of the slice.
%
%   Outputs:
%       shifts: A 3xN matrix where each column is the [x; y; 0] coordinate of the bottom-left corner of a cell for specified slice.

    offsets = double(slice) - 1;

    % Compute the slice cell starting coordinates.
    xIndices = mod(offsets, double(countX));
    yIndices = floor(double(offsets) / double(countX));

    xCoords = xIndices * sizeX;
    yCoords = yIndices * sizeY;

    shifts = [xCoords; yCoords; zeros(1, numel(offsets))];
end
