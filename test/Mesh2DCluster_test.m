clc; clearvars;

currentFolder = fileparts(mfilename('fullpath'));
addpath(currentFolder, '..', 'src/');

countX = 3;
countY = 3;
sizeX = 0.1;
sizeY = 0.2;
origin = [0; 0; 0];
roll = pi / 4;


myRectMeshPrimitive = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive, origin, pi/3, pi/2, roll);
myRectMesh2 = mom.MeshRect(myRectMeshPrimitive, [1, 1, 1], 0, 0, 0);

myMesh1 = myRectMesh1.getAsMesh();
myMesh2 = myRectMesh2.getAsMesh();

myMeshes = mom.Mesh2DCluster();
myMeshes.loadMesh(myMesh1);
myMeshes.loadMesh(myMesh2);

mesh = myMeshes.getMesh(2);
disp(mesh.faces);
disp(myMeshes.meshRanges)