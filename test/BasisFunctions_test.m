clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 2;
countY = 2;
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];

tic;

% single Mesh2D
myRectMeshPrimitive = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh = mom.MeshRect(myRectMeshPrimitive, origin, 0, 0, 0);

myMesh = myRectMesh.getAsMesh();
basisFunc = mom.BasisFunctions(myMesh);

order = 2;

quadPoints = mom.QuadPoints(order);

basisFunc.generatePoints(myMesh, quadPoints);

basisFunc.getFaces()

qP = basisFunc.getQuadPoints();
rP = basisFunc.getRooftopPoints();

toc;