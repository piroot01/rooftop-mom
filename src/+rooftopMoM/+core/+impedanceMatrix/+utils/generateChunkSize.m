function [rows, cols] = generateChunkSize(elementCount, order)
% generateChunkSize: Computes a nearly square chunk size in blocks.
%
%   Inputs:
%       elementCount: Total number of finite elements (positive integer).
%       order: Quadrature order (positive integer), block size = order^4.
%
%   Outputs:
%       rows, cols: Dimensions of the chunk in blocks.

    % Compute number of elements per block and total blocks per chunk
    elementsPerBlock = order^4;
    totalBlocks = floor(double(elementCount) / double(elementsPerBlock));

    % Handle trivial cases.
    if totalBlocks <= 1
        rows = totalBlocks;
        cols = 1;
        return;
    end

    % Start from sqrt(totalBlocks) and search downward for a divisor.
    maxRow = floor(sqrt(totalBlocks));
    for r = maxRow:-1:1
        if mod(totalBlocks, r) == 0
            rows = r;
            cols = totalBlocks / r;
            return;
        end
    end
end
