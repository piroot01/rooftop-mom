clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 2;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

[u, v] = rooftopMoM.mesh.utils.generateRectLocalCoords(mesh.vertices{1}, mesh.faces{1})

localVertices = rooftopMoM.mesh.utils.generateRectLocalVertices(mesh.vertices{1}, u, v)
