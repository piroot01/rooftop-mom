function chunks = generateChunks(blockCount, count)
% generateChunks Partition 1:count into contiguous chunks of size up to blockCount.
%
%   Inputs
%       blockCount: positive integer scalar. Maximum number of indices per chunk (i.e. block size).
%       count: nonnegative integer scalar. Total number of indices to partition (the range 1:count).
%
%   Output:
%       chunks: 1 x chunkCount cell array. Each cell is a struct with the following fields:
%           rowIndices: vector of consecutive indices for this chunk (startIdx:endIdx).
%           columnIndices: scalar (currently set to  value 1 by this implementation).

    chunkCount = ceil(double(count) / blockCount);
    chunks = cell(1, chunkCount);

    for n = 1:chunkCount
        startIdx = (n - 1) * blockCount + 1;
        endIdx = min(n * blockCount, count);

        chunks{n} = struct( ...
            'rowIndices', startIdx:endIdx, ...
            'columnIndices', 1 ...
        );
    end
end
