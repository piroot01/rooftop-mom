function area = generateRectFaceArea(vertices, faces, faceIndex)
% generateRectFaceArea: Returns an area of rectangular face.
%
%   Inputs:
%       vertices: 3×N array, each column is a 3D vertex coordinate.
%       faces: 4×F integer array, each column lists vertex indices of a quad face.
%       faceIndex: scalar index (indices) of the face(s).
%
%   Output:
%       origin: 3xN array of the origins of the selected face(s).

    % Default to the first face if none specified.
    if nargin < 3 || isempty(faceIndex)
        faceIndex = 1;
    end

    % Retrieve principal edges.
    [e1, e2] = rooftopMoM.mesh.utils.generateRectFacePrincipalEdges(vertices, faces, faceIndex);

    area = norm(cross(e1, e2));
end
