function [bMin, bMax] = generateRectBoundary(vertices, u, v)
% generateRectBoundary: Compute axis-aligned bounding rectangle in local face basis.
%
%   Inputs:
%       vertices: 3×N array, each column is a 3D vertex coordinate.
%       u, v: 3×1 vectors that serve as the local projection axes.
%
%   Output:
%       bMin, bMax: two 3d vectors, first one is the minimum corner and the second one is the maximum corner of the face's boundary expressed in global coordinates.

    % Project all vertices onto the local axes.
    localVertices = rooftopMoM.mesh.utils.generateRectLocalVertices(vertices, u, v);

    us = localVertices(1, :);
    vs = localVertices(2, :);

    % Determine min and max along each axis.
    uMin = min(us);
    uMax = max(us);
    vMin = min(vs);
    vMax = max(vs);

    % Map boundary corners back to global coordinates.
    bMin = uMin * u + vMin * v;
    bMax = uMax * u + vMax * v;
end
