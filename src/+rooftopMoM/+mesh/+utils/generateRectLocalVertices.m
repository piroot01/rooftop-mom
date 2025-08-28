function localVertices = generateRectLocalVertices(vertices, u, v)
% generateRectLocalVertices: Project 3D points onto a planar mesh.
%
%   Inputs:
%       vertices: 3×N array; each column is a 3D vertex coordinate.
%       u, v: 3×1 vectors that serve as the local projection axes.
%
%   Output:
%       localVertices: 2×N array; each column is the corresponding 2D projected vertex.

    projectors = [u(:), v(:)];
    localVertices = rooftopMoM.mesh.utils.generateProjectedVertices(vertices, projectors);
end
