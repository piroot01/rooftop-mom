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

face1Idx = mesh.faces(:, 1);
face1Vert = mesh.vertices(:, face1Idx)

v1 = face1Vert(:, 1)
v2 = face1Vert(:, 2)
v4 = face1Vert(:, 4)

order = 3;

myQuadPoints = mom.QuadPoints(order);
myQuadPoints.getTransformed(v1, v2, v4)
