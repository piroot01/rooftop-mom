function projectedVertices = generateProjectedVertices(vertices, projectors)
% generateRectProjectedVertices: Project 3D points onto multiple axes.
%
%   Inputs:
%       vertices: 3×N array, each column is a 3D point.
%       projectors: 3×M array, each column is a 3D projector (axis).
%
%   Output:
%       projectedVertices: M×N array: 
%           row i is the projection of all points onto projector i
%
%   Example:

    projectedVertices = projectors' * vertices;
end
