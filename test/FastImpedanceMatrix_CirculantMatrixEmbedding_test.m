clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 2;
countY = 2;
sizeX = 1;
sizeY = 1;
gridSize = [countX, countY];
origin = [0, 0, 0];

order = 5;

lambda = 5;
k = 2 * pi / lambda;

myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin, 0, 0, 0);
myMeshCluster = mom.Mesh2DCluster();
myMeshCluster.loadMesh(myRectMesh1.getAsMesh());

quadrature = mom.Quadrature(order);

basisFunctions = mom.BasisFunctions(myMeshCluster);

rooftopBasisFunctions = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize);
rooftopBasisFunctions.generatePoints();

strategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, rooftopBasisFunctions, k);
strategyMag = mom.FastImpedanceMatrix.ComputeMagPart();
strategyEl = mom.FastImpedanceMatrix.ComputeElPart();


tic;

[GT0x, ~] = strategyMag.compute(strategyCtx);
GP0 = strategyEl.compute(strategyCtx);

toc;

tic;

GT0xP = mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(GT0x, gridSize)
GP0P = mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(GP0, gridSize)

toc;
