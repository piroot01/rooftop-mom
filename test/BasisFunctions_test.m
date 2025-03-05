clc; clearvars;

currentFolder = fileparts(mfilename('fullpath'));
addpath(currentFolder, '..', 'src/');

countX = 10;
countY = 3;
sizeX = 0.1;
sizeY = 0.2;
origin = [0, 0, 0];

% single Mesh2D
myRectMeshPrimitive = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh = mom.MeshRect(myRectMeshPrimitive, origin, 0, 0, 0);

myMesh = myRectMesh.getAsMesh();
basisFunc = mom.BasisFunctions(myMesh);
basisFunc.get()

% multiple mesh2D -> mesh2DCluster
myMeshCluster = mom.Mesh2DCluster();
myRectMeshShifted = mom.MeshRect(myRectMeshPrimitive, origin + [0, 0, 10], 0, 0, 0);
myMeshCluster.loadMesh(myRectMeshShifted.getAsMesh());
myMeshCluster.loadMesh(myMesh);

basisFunc = mom.BasisFunctions(myMeshCluster);
basisFunc.get()
