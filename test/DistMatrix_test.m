clc; clearvars;

currentFolder = fileparts(mfilename('fullpath'));
addpath(currentFolder, '..', 'src/');

countX = 3;
countY = 3;
sizeX = 0.5;
sizeY = 1.5;
origin = [0, 0, 2];

rectPrim = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
rectMesh = mom.MeshRect(rectPrim, origin, 0, 0, 0);

mesh = rectMesh.getAsMesh();

vertexIndices = mesh.faces(:, [1, 2, 3, 4]);
faces = mesh.vertices(:, vertexIndices(:))

order = 2;

myQuadPoints = mom.QuadPoints(order);
quadPoints = myQuadPoints.getTransformed(faces)

mom.DistMatrix().getChunk(quadPoints(:, 1:(4 * 2)), quadPoints(:, 9:(4 * 4)))

A = magic(5);
A(2:3, 2:4) = 1:6