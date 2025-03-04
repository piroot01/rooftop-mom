clc; clearvars;

currentFolder = fileparts(mfilename('fullpath'));
addpath(currentFolder, '..', 'src/');

countX = 10;
countY = 3;
sizeX = 0.1;
sizeY = 0.2;
origin = [0; 0; 0];
roll = pi / 4;


myRectMeshPrimitive = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh = mom.MeshRect(myRectMeshPrimitive, origin, pi/3, pi/2, roll);

myMesh = myRectMesh.getAsMesh();
basisFunc = mom.BasisFunctions(myMesh);
basisFunc.get()