function memory = memoryForRooftopPoints(rooftops)
% memoryForRooftopPoints: Returns the memory footprint of rooftop point data.
%
%   INPUT:
%       rooftops: Instance of rooftopMoM.basisFunctions.(Mesh/UniformRectPrimitive)UniformRooftops.
%
%   OUTPUT:
%       memory: Total memory in bytes needed for storing both positive and negative rooftopPoints and quadraturePoints arrays.
%
%   NOTE: The memory does not include the cell memory overhead.

    % Compute number of 3D points per rooftop: order^2 nodes times 3 coordinates.
    singleRooftopCount = rooftops.quadrature().order()^2 * 3;

    % Total number of rooftop elements.
    count = rooftops.count();

    % There are 4 arrays: rooftopPointsPositive, rooftopPointsNegative, quadraturePointsPositive, quadraturePointsNegative.
    numArrays = 4;

    % Memory per double in bytes.
    bytesPerDouble = rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.Double);

    % Total memory footprint.
    memory = numArrays * count * singleRooftopCount * bytesPerDouble;
end
