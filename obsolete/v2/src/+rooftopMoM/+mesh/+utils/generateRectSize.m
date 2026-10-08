function [s1, s2] = generateRectSize(vertices, u, v)
% generateRectGridSize: Determine geometrical sizes along principal directions.
%
%   Inputs:
%       vertices: 3×N array, each column is a 3D vertex coordinate.
%       u, v: 3×1 vectors that serve as the local projection axes.
%
%   Outputs:
%       s1, s2: Sizes of the rectangular patch boundaries along each principal edges.

    % Project all vertices onto the local axes.
    localVertices = rooftopMoM.mesh.utils.generateRectLocalVertices(vertices, u, v);

    us = localVertices(1, :);
    vs = localVertices(2, :);

    % Determine extents along each axis.
    uMin = min(us);
    uMax = max(us);
    vMin = min(vs);
    vMax = max(vs);

    s1 = uMax - uMin;
    s2 = vMax - vMin;
end
