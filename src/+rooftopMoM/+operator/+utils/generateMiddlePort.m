function port = generateMiddlePort(rooftops, direction)
% generateMiddlePort: Returns the port across the middle of the grid, the rooftops of one direction on the middle
%   line of edges across it.
%
%   Inputs:
%       rooftops: rooftopMoM.basisFunctions.UniformRooftops with rooftops of the direction.
%       direction: 'x' or 'y', the direction of the current through the port.
%
%   Output:
%       port: 1xP indices of the rooftops of the port.
%
%   NOTE: An odd number of faces has no line of edges in the middle; the line just before it is taken.

    primitive = rooftops.primitive();
    [rowCount, columnCount] = rooftops.gridSize(direction);

    % The middle line is one row of the x-directed rooftops or one column of the y-directed ones.
    if direction == 'x'
        rows = repmat(floor(primitive.countX() / 2), 1, columnCount);
        columns = 1:columnCount;
        first = 0;
    else
        rows = 1:rowCount;
        columns = repmat(floor(primitive.countY() / 2), 1, rowCount);
        first = rooftops.countDir('x');
    end

    port = first + sub2ind([rowCount, columnCount], rows, columns);
end
