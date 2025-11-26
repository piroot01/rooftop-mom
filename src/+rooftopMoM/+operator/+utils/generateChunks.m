function [chunks, chunkCount] = generateChunks(vectorSize, blockSize, chunkSize)
    % Single chunk determines an array of size 3 x blockSize
    blockCountPerOneChunk = max(1, ceil(double(chunkSize) / double(blockSize)));

    % Compute start and end indices for each chunk.
    starts = 1:blockCountPerOneChunk:vectorSize;
    ends = starts + blockCountPerOneChunk - 1;
    ends(end) = min(ends(end), vectorSize);

    chunkCount = numel(starts);
    chunks = cell(1, chunkCount);

    for n = 1:chunkCount
        chunks{n} = struct( ...
            'indicesStart', starts(n), ...
            'indicesEnd',   ends(n) ...
        );
    end
end
