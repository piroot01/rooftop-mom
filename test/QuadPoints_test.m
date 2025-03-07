clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 0.5;
sizeY = 1.5;
origin = [0, 0, 2];

rectPrim = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
rectMesh = mom.MeshRect(rectPrim, origin, 0, 0, 0);

mesh = rectMesh.getAsMesh();

vertexIndices = mesh.faces(:, [1, 2, 3]);
faces = mesh.vertices(:, vertexIndices(:))

order = 2;

myQuadPoints = mom.QuadPoints(order);
myQuadPoints.getTransformed(faces)