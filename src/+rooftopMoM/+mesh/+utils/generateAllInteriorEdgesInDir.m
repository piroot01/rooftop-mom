function indices = generateAllInteriorEdgesInDir(countX, countY, direction)
% generateAllInteriorEdgesInDir: Return face-index pairs for grid edges in a given direction.
%
%   Inputs:
%       countX: Number of cells in the x-direction.
%       countY: Number of cells in the y-direction.
%       direction: 'x' or 'y' (string) specifying which edge direction to produce.
%
%   Output:
%       indices: 2-by-N array. Each column is a pair of linear face indices
%           [faceA; faceB] that share an edge. The pair ordering follows the adjacency along the requested grid dimension.
%
%   Indexing convention:
%       Linear face indices are arranged columnwise into a cX-by-cY matrix.

    % Build grid of linear face indices.
    idx = reshape(1:(countX * countY), countX, countY);

    switch direction
        case 'x'
            % Face pairs adjacent along the first grid dimension (x).
            % For each column, pair index (r, c) with (r+1, c).
            leftIdx  = idx(1:(end - 1), :);
            rightIdx = idx(2:end, :);
            indices = [leftIdx(:).'; rightIdx(:).'];
        case 'y'
            % Face pairs adjacent along the second grid dimension (y).
            % For each row, pair index (r, c) with (r, c+1).
            bottomIdx = idx(:, 1:(end - 1));
            topIdx    = idx(:, 2:end);
            indices = [bottomIdx(:).'; topIdx(:).'];
        otherwise
            error("Unknown direction '%s'. Use 'x' or 'y'.", direction);
    end
end
