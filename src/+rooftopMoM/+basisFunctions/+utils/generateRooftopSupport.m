function support = generateRooftopSupport(vertices, faces, facePair)
% generateRooftopMesh: Returns array of vertices forming a single rooftop defined of a pair of quads.
%
%   Inputs:
%       vertices: 3×N array, each column is a 3D vertex coordinate.
%       faces: 4×F integer array, each column lists vertex indices of a quad face.
%       facePair: Pair of two indices that forms the rooftop.
%
%   Output:
%       support: 3x8 vertex array fully defining the support of the rooftop.

    vertexIndices = faces(:, facePair);
    support = vertices(:, vertexIndices(:));
end
