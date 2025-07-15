function indices = generateSubmatrixIndices(indices, includeDiagonal)
% generateSubmatrixIndices: Returns an array of submatrix indices to include in future computations.
%
%   Inputs:
%       indices: 1×M cell array; each cell is a 2×K array of basis-function indices
%       includeDiagonal: Logical flag that determines whether diagonal
%                        submatrices should be included.
%
%   Output:
%       indices: 1×N array of submatrix indices to include.

    count = numel(indices);
    submatrixCount = count^2;

    if includeDiagonal || count == 1
        indices = 1:submatrixCount;
    else
        allIndices = 1:submatrixCount;
        diagIndices = 1:(count + 1):submatrixCount;
        indices = setdiff(allIndices, diagIndices);
    end
end
