function count = generateRectFaceCount(vertices, faces, faceIndex)
% generateRectFaceCount: Returns an number of faces in the full rect.
%
%   Inputs:
%       vertices: 3×N array, each column is a 3D vertex coordinate.
%       faces: 4×F integer array, each column lists vertex indices of a quad face.
%       faceIndex: scalar index (indices) of the face(s).
%
%   Output:
%       count: Number of the faces in the full rect.

    % Default to the first face if none specified.
    if nargin < 3 || isempty(faceIndex)
        faceIndex = 1;
    end

    [cX, cY] = rooftopMoM.mesh.utils.generateRectGridSize(vertices, faces, faceIndex);

    count = times(cX, cY);
end
