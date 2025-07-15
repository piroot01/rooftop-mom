function [rows, cols] = generateChunkSize(elementCount, order)
% generateChunkSize: Computes a nearly square chunk size in blocks.
%
%   Inputs:
%       elementCount: Total number of finite elements (positive integer).
%       order: Quadrature order (positive integer), block size = order^4.
%
%   Outputs:
%       rows, cols: Dimensions of the chunk in blocks.

    blockCount = rooftopMoM.core.impedanceMatrix.utils.generateBlockCountPerChunk(elementCount, order);

    % Handle trivial cases.
    if blockCount <= 1
        rows = blockCount;
        cols = 1;
        return;
    end

    % Start from sqrt(blockCount) and search downward for a divisor.
    maxRow = floor(sqrt(blockCount));
    for r = maxRow:-1:1
        if mod(blockCount, r) == 0
            rows = r;
            cols = blockCount / r;
            return;
        end
    end
end
