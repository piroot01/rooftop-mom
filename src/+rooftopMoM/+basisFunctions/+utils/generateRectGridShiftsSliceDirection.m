function shifts = generateRectGridShiftsSliceDirection(countX, sizeX, sizeY, slice, direction)
% generateRectGridShifts: Generates shift vectors for a rectangular grid.
%
%   Inputs:
%       countX: Number of cells in the x-direction.
%       countY: Number of cells in the y-direction.
%       sizeX:  Width of each cell in the x-direction.
%       sizeY:  Height of each cell in the y-direction.
%       slice: Indices of the slice.
%       direction: Direction of the shift indexing, 'x' or 'y'.
%
%   Outputs:
%       shifts: A 3xN matrix where each column is the [x; y; 0] coordinate of the bottom-left corner of a cell for specified slice.

    switch direction
        case 'x'
            shifts = rooftopMoM.basisFunctions.utils.generateRectGridShiftsSlice( ...
                countX - 1, sizeX, sizeY, ...
                slice ...
            );
        case 'y'
            shifts = rooftopMoM.basisFunctions.utils.generateRectGridShiftsSlice( ...
                countX, sizeX, sizeY, ...
                slice ...
            );
        otherwise
            error("Unknown direction '%s'. Use 'x' or 'y'.", direction);
    end
end
