function [offsets, sizes] = generateSubmatrixOffsetsAndSizes(indices)
% generateSubmatrixOffsetsAndSizes: Precompute submatrix offsets and sizes.
%
%   Input:
%       indices: 1×M cell array; each cell is a 2×K array of basis-function indices
%           (rows correspond to polarities 1 and 2; columns correspond to positions).
%
%   Outputs:
%       offsets: 2×N array of starting (row, column) offsets in the impedance matrix for each submatrix.
%       sizes: 2×N array of submatrix sizes (number of rows, number of columns).

    counts = cellfun(@(A) size(A, 2), indices);
    count = numel(counts);

    % Compute starting offsets along the diagonal for each block.
    diagonalOffsets = [0, cumsum(counts(1:(end - 1)))];

    % Generate all combinations of block indices.
    [I, J] = ndgrid(1:count);

    % Map block indices to matrix offsets and sizes.
    offsets = [diagonalOffsets(I(:)); diagonalOffsets(J(:))];
    sizes   = [counts(I(:)); counts(J(:))];
end
