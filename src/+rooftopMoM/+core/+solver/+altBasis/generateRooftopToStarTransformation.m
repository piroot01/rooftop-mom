function matrix = generateRooftopToStarTransformation(starsConverted, sizeX, sizeY)
% generateRooftopToStarTransformation: Generate a sparse matrix describing the transformation
%   from the rooftop basis into the star basis.
%
%   Inputs:
%       starsConverted: 4xN matrix containing rooftop indices for each patch/star.
%       sizeX: Size in the x-dir of a single rectangular patch.
%       sizeY: Size in the y-dir of a single rectangular patch.
%
%   Output:
%       matrix: Sparse matrix representation of the transformation.

    % The number of star basis functions is simply the number of patches (columns)
    numStars = size(starsConverted, 2);

    % The total number of rooftops.
    numRooftops = starsConverted(4, end);

    % Duplicate row coordinates for all 4 edges (Right, Top, Left, Bottom)
    rowIdx = repmat(1:numStars, 4, 1);

    % Column indices are the rooftop numbers directly from the input matrix
    colIdx = starsConverted(1:4, :);

    % Expand the values to match the size of our coordinates.
    % NOTE: Star basis functions typically represent divergence (flow out of the patch).
    % Adjust these signs/scalings if your specific rooftop orientation convention 
    % (e.g., inward pointing vs outward pointing normal vectors) dictates otherwise.
    vals = repmat([1/sizeX; 1/sizeY; -1/sizeX; -1/sizeY], 1, numStars);

    % Filter out zero indices (boundary edges with no neighboring patch)
    valid = colIdx > 0;

    % Assemble the sparse matrix in one shot
    matrix = sparse(rowIdx(valid), colIdx(valid), vals(valid), numStars, numRooftops);
end
