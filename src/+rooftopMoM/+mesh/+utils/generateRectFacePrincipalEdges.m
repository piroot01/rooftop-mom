function [e1, e2] = generateRectFacePrincipalEdges(vertices, faces, faceIndex)
% extractRectEdges: Retrieve the two principal edge vectors of a quad face.
%
%   Inputs:
%       vertices: 3×N array, each column is a 3D vertex coordinate.
%       faces: 4×F integer array, each column lists vertex indices of a quad face.
%       faceIndex: (optional) scalar index of the face (default: 1).
%
%   Output:
%       e1, e2: 3×1 vectors for edges from the first vertex to the second (e1) and to the fourth (e2) of the specified face.

    % Default to the first face if none specified.
    if nargin < 3 || isempty(faceIndex)
        faceIndex = 1;
    end

    % Extract indices of the selected face.
    idx = faces(:, faceIndex);

    % Fetch corner vertices: v1, v2, v4.
    v1 = vertices(:, idx(1));
    v2 = vertices(:, idx(2));
    v4 = vertices(:, idx(4));

    % Compute and return edge vectors.
    e1 = v2 - v1;
    e2 = v4 - v1;
end
