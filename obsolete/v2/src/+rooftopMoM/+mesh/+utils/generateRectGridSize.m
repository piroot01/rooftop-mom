function [c1, c2] = generateRectGridSize(vertices, faces, faceIndex)
% generateRectGridSize: Determine grid counts along principal directions.
%
%   Input:
%       vertices: 3×N array, each column is a 3D vertex coordinate.
%       faces: 4×F integer array of quad face vertex indices.
%       faceIndex: (optional) scalar index of the face (default: 1).
%
%   Output:
%       c1, c2: Unsigned integer counts of subdivisions along the u- and v-axes, computed by dividing the face extents by the respective edge lengths.

    % Default to the first face if unspecified.
    if nargin < 3 || isempty(faceIndex)
        faceIndex = 1;
    end

    [u, v] = rooftopMoM.mesh.utils.generateRectLocalCoords(vertices, faces);

    % Compute plate boundary sizes.
    [s1, s2] = rooftopMoM.mesh.utils.generateRectSize(vertices, u, v);

    % Get face edge lengths.
    [l1, l2] = rooftopMoM.mesh.utils.generateRectFaceSize(vertices, faces, faceIndex);

    % Compute subdivision counts and convert to uint32.
    c1 = uint32(round(s1 / l1));
    c2 = uint32(round(s2 / l2));
end
