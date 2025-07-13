function origin = generateRectFaceOrigin(vertices, faces, faceIndex)
% generateRectFaceOrigin: Obtain the origin(s) of a selected quad face(s).
%
%   Inputs:
%       vertices: 3×N array, each column is a 3D vertex coordinate.
%       faces: 4×F integer array, each column lists vertex indices of a quad face.
%       faceIndex: scalar index (indices) of the face(s).
%
%   Output:
%       origin: 3xN array of the origins of the selected face(s).

    % Extract indices of the selected face.
    idx = faces(:, faceIndex(:));

    origin = vertices(:, idx(1, :));
end
