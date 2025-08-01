function transformed = transformIndices(indices, fullSize, blockSize, rowOffset, columnOffset, extraShift)
% transformIndices: Maps local block indices to global matrix indices.
%   Converts a list of linear indices ('blkIndices') that are local to a square block into linear indices valid within a larger global matrix.
%
%   The transformation applies both a geometric placement (via row/column
%   offsets) and a final numerical shift.
%
%   Inputs:
%       indices:  A vector of linear indices defined within the local block.
%       fullSize: A 1x2 vector [rows, cols] specifying the global matrix dimensions.
%       blockSize: A 1x2 vector [rows, cols] dimension of the local square block.
%       rowOffset: The row offset to apply for geometric placement (moves block down).
%       columnOffset: The column offset to apply for geometric placement (moves block right).
%       extraShift: A scalar value added to the final global indices after placement.
%
%   Output:
%       transformed: A vector of the resulting global linear indices, cast as uint32.

    if nargin < 6 || isempty(extraShift)
        extraShift = 0;
    end

    % Convert local linear indices to local [row, col] subscripts.
    [r, c] = ind2sub(blockSize, indices);

    % Apply geometric offsets to translate local subscripts to global positions.
    r = r + rowOffset;
    c = c + columnOffset;

    % Convert global subscripts back to linear indices in the full-size matrix,
    transformed = uint32(sub2ind(fullSize, r, c));

    % Shift the indices by extraShift.
    transformed = rooftopMoM.core.impedanceMatrix.utils.shiftIndices(transformed, extraShift);
end
