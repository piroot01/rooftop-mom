function slices = generateSlicesFromChunks(chunks, order)
% generateSlicesFromChunks: Converts block indices (chunks) to linear indices (slices).
%
%   Inputs:
%       chunks: An Mx1 cell array. Each cell must contain a struct with the fields 'rowIndices' and 'columnIndices'
%       order:  A scalar integer representing the quadrature order. This determines the size of each block.
%
%   Output:
%       slices: An Mx1 cell array. Each cell contains a struct defining the start and end linear indices for rows and columns. 
%               The struct fields are 'rowIndicesStart', 'rowIndicesEnd', 'columnIndicesStart', and 'columnIndicesEnd'.

    % The size of each fine-grained block is the square of the quadrature order.
    blockSize = order^2;

    slices = cellfun(@(c) struct(...
        'rowIndicesStart', (c.rowIndices(1) - 1) * blockSize + 1, ...
        'rowIndicesEnd', c.rowIndices(end) * blockSize, ...
        'columnIndicesStart', (c.columnIndices(1) - 1) * blockSize + 1, ...
        'columnIndicesEnd', c.columnIndices(end) * blockSize), ...
        chunks, 'UniformOutput', false);
end
