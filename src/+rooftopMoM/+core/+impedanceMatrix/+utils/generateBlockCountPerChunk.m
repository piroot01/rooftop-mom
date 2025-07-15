function count = generateBlockCountPerChunk(elementCount, order)
% generateBlockCountPerChunk: Determines the number of blocks (elementary units) per chunk.
%
%   Inputs:
%       elementCount: Total number of finite elements (positive integer).
%       order:        Quadrature order (positive integer); block size = order^4.
%
%   Output:
%       count:        Number of blocks per chunk.

    elementsPerBlock = order^4;

    if elementCount < elementsPerBlock
        error('Need at least %d elements worth of memory, only %d available.', elementsPerBlock, elementCount);
    end

    count = floor(double(elementCount) / double(elementsPerBlock));
end
