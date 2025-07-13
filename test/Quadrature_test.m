clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 2;
countY = 2;
sizeX = 0.5;
sizeY = 1.5;
origin = [0, 0, 2];

rectPrim = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
rectMesh = mom.MeshRect(rectPrim, origin, 0, 0, 0);

mesh = rectMesh.getAsMesh();

vertexIndices = mesh.faces(:, 1:2);
vertexCoords = mesh.vertices(:, vertexIndices(:));

order = 2;
myQuadrature = mom.Quadrature(order);
myQuadrature.get(vertexCoords)
