clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 2;
countY = 2;
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];

% single Mesh2D
myRectMeshPrimitive = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh = mom.MeshRect(myRectMeshPrimitive, origin, 0, 0, 0);

myMesh = myRectMesh.getAsMesh();
basisFunc = mom.BasisFunctions(myMesh);

order = 2;

quadPoints = mom.QuadPoints(order);

basisFunc.generatePoints(myMesh, quadPoints);

memoryElementwise = 16;

memory_MB = memoryElementwise * 8 / 1e6;

impedMatrix = mom.ImpedanceMatrixChunked(memory_MB, order, basisFunc);

impedMatrix.compute(basisFunc, quadPoints.getWeights4D());