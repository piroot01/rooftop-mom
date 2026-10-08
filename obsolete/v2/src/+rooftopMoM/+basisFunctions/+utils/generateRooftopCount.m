function count = generateRooftopCount(countX, countY)
% generateRooftopCount: Return the maximum number of rooftops on a grid.
%
%   Inputs:
%       countX, countY: Grid dimensions for the antenna topology.
%
%   Output:
%       count: Maximum number of rooftops that can occupy the grid.

    [cX, cY] = rooftopMoM.basisFunctions.utils.generateRooftopCounts(countX, countY);

    count = cX + cY;
end
