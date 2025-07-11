function indices = filterInteriorEdgesByDir(interiorEdgesIndices, dir)
% filterInteriorEdgesByDir: Filter interior mesh‐edge face pairs by the edge axis direction.
%
%   Inputs:
%       interiorEdgesIndices: 2×E array, where each column is a pair of face indices sharing an interior edge.
%       dir: Character scalar, either:
%           'x' – select edges where face indices differ by exactly 1,
%           'y' – select all other edges.
%
%   Output:
%       indices: 2×E_d array of face‐index pairs for which the edge orientation matches the requested axis.

    f1 = interiorEdgesIndices(1, :);
    f2 = interiorEdgesIndices(2, :);

    % Identify edges aligned in the 'x' direction (adjacent face indices differ by 1).
    diffs = abs(int64(f1) - int64(f2));
    isXDir = (diffs == 1);

    % Return edges matching the requested direction.
    if dir == 'x'
        indices = interiorEdgesIndices(:, isXDir);
    else
        indices = interiorEdgesIndices(:, ~isXDir);
    end
end
