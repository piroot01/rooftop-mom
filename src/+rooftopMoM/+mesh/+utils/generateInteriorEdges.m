function indices = generateInteriorEdges(vertices, faces, tol)
% generateInnerEdges: Finds all interior edges of a quad mesh.
%
%   input:
%       faces: 4xN array of vertex indices (each comlumn is one quad)
%   output:
%       indices: 2xE array of face-pairs sharing the interior edge
    if nargin < 3
        tol = eps;
    end

    faceCount = size(faces, 2);

    % Build a local 2D frame in the mesh plane.
    v1 = vertices(:, faces(1, 1))
    v2 = vertices(:, faces(2, 1))
    v4 = vertices(:, faces(4, 1))
    e1 = v2 - v1;
    e2 = v4 - v1;

    % Create an orthonormal basis (u, v).
    u = e1 / norm(e1);

    % Create a rejection of the e2 to the span(u).
    e2Rejection = e2 - (dot(e2, u) * u);
    v = e2Rejection / norm(e2Rejection);

    % Project all vertices into the (u, v) coordinates.
    localVertices = vertices - v1;
    uCoords = u.' * vertices;
    vCoords = v.' * vertices;

    % Cluster unique grid-lines with tol tolerance.
    [uLines, ~, uIndices] = uniquetol(uCoords.', tol, 'ByRows', false)
    [vLines, ~, vIndices] = uniquetol(vCoords.', tol, 'ByRows', false)

    uLinesCount = numel(uLines)
    vLinesCount = numel(vLines)

    faceGrid = reshape(1:faceCount, uLinesCount - 1, vLinesCount - 1)
end
