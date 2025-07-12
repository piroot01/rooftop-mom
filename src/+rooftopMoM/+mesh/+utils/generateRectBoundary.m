function [bMin, bMax] = generateRectBoundary(vertices, faces)
% generateRectBoundary: Compute axis-aligned bounding rectangle in local face basis.
%
%   Inputs:
%       vertices: 3×N array, each column is a 3D vertex coordinate.
%       faces: 4×F integer array, each column lists vertex indices of a quad face.
%
%   Output:
%       bMin, bMax: two 3d vectors, first one is the minimum corner and the second one is the maximum corner of the face's boundary expressed in global coordinates.

    % Compute local basis for the face.
    [u, v] = rooftopMoM.mesh.utils.generateRectLocalCoords(vertices, faces);

    % Project all vertices onto the local axes.
    us = u.' * vertices;
    vs = v.' * vertices;

    % Determine min and max along each axis.
    uMin = min(us);
    uMax = max(us);
    vMin = min(vs);
    vMax = max(vs);

    % Map boundary corners back to global coordinates.
    bMin = uMin * u + vMin * v;
    bMax = uMax * u + vMax * v;
end
