function indices = generateInteriorEdges(faces)
% generateInnerEdges: Finds all interior edges of a quad mesh.
%
%   Input:
%       faces: 4xN array of vertex indices (each comlumn is one quad)
%   Output:
%       indices: 2xE array of face-pairs sharing the interior edge

    faceCount = size(faces, 2);

    v1 = faces(1, :).';
    v2 = faces(2, :).';
    v3 = faces(3, :).';
    v4 = faces(4, :).';

    % Build all edges with sorted vertices.
    faceIndices   = uint32(1:faceCount).';
    edges = [ min(v1, v2), max(v1, v2), faceIndices; ...
              min(v2, v3), max(v2, v3), faceIndices; ...
              min(v3, v4), max(v3, v4), faceIndices; ...
              min(v4, v1), max(v4, v1), faceIndices ];

    % Group edges by vertex pairs.
    [~, ~, ic] = unique(edges(:, 1:2), 'rows');

    % Count occurrences per edge.
    counts = accumarray(ic, 1);

    % Collect face indices for each edge and sort them.
    facePairs = accumarray(ic, edges(:, 3), [], @(xs) {sort(xs)});

    % Select edges shared by exactly two faces (interior).
    interiorEdgesMask = (counts == 2);
    interiorFacePairs = facePairs(interiorEdgesMask);

    % Convert to 2xE output.
    indices = [interiorFacePairs{:}];
end
