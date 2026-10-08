function [l1, l2] = generateRectFaceSize(vertices, faces, faceIndex)
% generateRectFaceSize: Compute the side lengths of a quadrilateral face.
%
%   Inputs:
%       vertices: 3×N array, each column is a 3D vertex coordinate.
%       faces: 4×F integer array, each column lists vertex indices of a quad face.
%       faceIndex: (optional) scalar index of the face (default: 1).
%
%   Output:
%       l1, l2: lengths of edges from the first vertex to the second and from the first vertex to the fourth of the specified face.

    % Default to the first face if none specified.
    if nargin < 3 || isempty(faceIndex)
        faceIndex = 1;
    end

    % Retrieve principal edges.
    [e1, e2] = rooftopMoM.mesh.utils.generateRectFacePrincipalEdges(vertices, faces, faceIndex);

    l1 = norm(e1);
    l2 = norm(e2);
end
