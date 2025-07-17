function [s1, s2] = generateRectSize(vertices, faces, faceIndex)
% generateRectGridSize: Determine geometrical sizes along principal directions.
%
%   Input:
%       vertices: 3×N array, each column is a 3D vertex coordinate.
%       faces: 4×F integer array of quad face vertex indices.
%       faceIndex: (optional) scalar index of the face (default: 1).
%
%   Output:
%       s1, s2: Sizes of the rectangular patch boundaries along each principal edges.

    % Default to the first face if unspecified.
    if nargin < 3 || isempty(faceIndex)
        faceIndex = 1;
    end

    % Compute local basis vectors.
    [u, v] = rooftopMoM.mesh.utils.generateRectLocalCoords(vertices, faces, faceIndex);

    % Project vertices onto local axes.
    us = u.' * vertices;
    vs = v.' * vertices;

    % Determine extents along each axis.
    uMin = min(us);
    uMax = max(us);
    vMin = min(vs);
    vMax = max(vs);

    s1 = uMax - uMin;
    s2 = vMax - vMin;
end
