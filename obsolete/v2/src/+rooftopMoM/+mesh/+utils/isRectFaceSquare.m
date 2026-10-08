function verdict = isRectFaceSquare(vertices, faces, faceIndex)
% isRectFaceSquare: Determine whether a rectangular face is square.
%
%   Inputs:
%       vertices: 3×N array, where each column is a 3D vertex coordinate.
%       faces: 4×F integer array, where each column lists the vertex indices of a quadrilateral face.
%       faceIndex: (optional) Index of the face to test (default: 1).
%
%   Output:
%       verdict: Logical value indicating whether the specified face is square.

    % Default to the first face if none specified.
    if nargin < 3 || isempty(faceIndex)
        faceIndex = 1;
    end

    [l1, l2] = rooftopMoM.mesh.utils.generateRectFaceSize(vertices, faces, faceIndex);

    verdict = (l1 == l2);
end
