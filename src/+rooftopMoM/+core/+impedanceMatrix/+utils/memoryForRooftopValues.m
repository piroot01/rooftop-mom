function memory = memoryForRooftopValues(rooftops)
% memoryForRooftopPoints: Returns the memory footprint of rooftop value data.
%
%   INPUT:
%       rooftops: Instance of rooftopMoM.basisFunctions.(Mesh/UniformRectPrimitive)UniformRooftops.
%
%   OUTPUT:
%       memory: Total memory in bytes needed for storing both positive and negative rooftopDivPoints and rooftop dimensions.
%
%   NOTE: The memory does not include the cell memory overhead.

    % Total number of rooftop elements.
    count = rooftops.count();

    % There are 4 arrays: rooftopDivPointsPositive, rooftopDivPointsNegative, length, height.
    numArrays = 4;

    % Memory per double in bytes.
    bytesPerDouble = rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.Double);

    % Total memory footprint.
    memory = numArrays * count * bytesPerDouble;
end

