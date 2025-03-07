clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 2;
countY = 2;
sizeX = 0.1;
sizeY = 0.2;
origin = [0, 0, 0];

% single Mesh2D
myRectMeshPrimitive = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh = mom.MeshRect(myRectMeshPrimitive, origin, 0, 0, 0);

myMesh = myRectMesh.getAsMesh();
basisFunc = mom.BasisFunctions(myMesh);

order = 2;

quadPoints = mom.QuadPoints(order);

basisFunc.generatePoints(myMesh, quadPoints)

[pointsPos, pointsNeg] = basisFunc.getPoints();

memoryElementwise = 48;

memory_MB = memoryElementwise * 8 / 1e6;

impedMatrix = mom.ImpedanceMatrixChunked(memory_MB, order, basisFunc);

impedMatrix.compute();