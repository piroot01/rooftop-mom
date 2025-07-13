function [u, v] = generateRectLocalCoords(vertices, faces, faceIndex)
% generateRectLocalCoords: Compute a local orthonormal basis for a given quadrilateral face.
%
%   Inputs:
%       vertices: 3×N array, columns are 3D vertex coordinates.
%       faces: 4×F integer array, columns are indices of vertices forming each quad face.
%       faceIndex: (optional) scalar index of the target face (default: 1).
%
%   Output:
%       u, v: ordered pair (u, v), two orthonormal basis vectors spanning the plane of the specified face.

    % Default to the first face if none specified.
    if nargin < 3 || isempty(faceIndex)
        faceIndex = 1;
    end

    % Retrieve principal edges.
    [e1, e2, ~] = rooftopMoM.mesh.utils.generateRectFaceFrame(vertices, faces, faceIndex);

    % First basis vector: normalized e1.
    u = e1 / norm(e1);

    % Second basis vector: component of e2 orthogonal to u.
    e2Proj = dot(e2, u) * u;
    e2Rej = e2 - e2Proj;
    v = e2Rej / norm(e2Rej);
end
