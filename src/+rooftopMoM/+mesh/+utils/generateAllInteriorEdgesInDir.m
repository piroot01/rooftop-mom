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

    switch direction
        case 'x'
            if countX < 2 || countY < 1
                indices = zeros(2,0);
                return
            end

            r = (1:(countX - 1)).';
            offs = (0:(countY - 1)) * countX;
            left = r + offs;
            left = left(:).';
            indices = cat(1, left, left + 1);
        case 'y'
            if countY < 2 || countX < 1
                indices = zeros(2,0);
                return
            end

            r = (1:countX).';
            offs = (0:(countY - 2)) * countX;
            bottom = r + offs;
            bottom = bottom(:).';
            indices = cat(1, bottom, bottom + countX);
        otherwise
            error("Unknown direction '%s'. Use 'x' or 'y'.", direction);
    end
end
