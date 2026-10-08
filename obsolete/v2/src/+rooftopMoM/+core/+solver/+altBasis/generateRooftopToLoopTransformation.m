function matrix = generateRooftopToLoopTransformation(loopsConverted, sizeX, sizeY)
% generateRooftopToLoopTransformation: Generate a sparse matrix describing the transformation
%   from the rooftop basis into the loop basis.
%
%   Inputs:
%       loopsConverted: Representation of the loops in the means of the rooftop basis.
%       sizeX: Size in the x-dir of the single rectangle.
%       sizeY: Size in the y-dir of the single rectangle.
%
%   Output:
%       matrix: Sparse matrix representation of the transformation matrix, filled with O(N)
%           elements.

    % The number of the rooftop basis functions will always be
    % contained inside the loopBasis matrix at the second row
    % last column
    numRooftops = loopsConverted(2, end);
    numLoops = size(loopsConverted, 2);

    % Duplicate row coordinates for all 4 edges
    rowIdx = repmat(1:numLoops, 4, 1);
    colIdx = loopsConverted(1:4, :);

    % Expand the values to match the size of our coordinates
    vals = repmat([1/sizeY; 1/sizeX; -1/sizeY; -1/sizeX], 1, numLoops);

    % Filter out invalid/zero indices (e.g., boundary edges)
    valid = colIdx > 0 & colIdx <= numRooftops;

    % Assemble the sparse matrix in one shot
    matrix = sparse(rowIdx(valid), colIdx(valid), vals(valid), numLoops, numRooftops);
end
