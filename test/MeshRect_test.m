clc; clearvars;

currentFolder = fileparts(mfilename('fullpath'));
addpath(currentFolder, '..', 'src/');

countX = 3;
countY = 3;
sizeX = 0.1;
sizeY = 0.2;
origin = mom.Vec3(0, 0, 0);
roll = pi / 4;


myRectMeshPrimitive = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh = mom.MeshRect(myRectMeshPrimitive, origin, pi/3, pi/2, roll);

myMesh = myRectMesh.getAsMesh();

disp(myMesh.vertices);
disp(myMesh.faces);