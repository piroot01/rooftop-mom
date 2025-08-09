function mesh = generateRectMesh(countX, countY, sizeX, sizeY)
% generateRectMesh: Creates and returns a single continuous mesh with no internal gaps or holes.
%
%   Inputs:
%       countX: Number of rectangular cells along the x-axis direction.
%       countY: Number of rectangular cells along the y-axis direction.
%       sizeX: Physical length of each cell along the x-axis.
%       sizeY: Physical length of each cell along the y-axis.
%
%   Output:
%       mesh: A rooftopMoM.mesh.Mesh object representing the generated uniform rectangular mesh.

    rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
    rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
    mesh = rect.generateMesh();
end
