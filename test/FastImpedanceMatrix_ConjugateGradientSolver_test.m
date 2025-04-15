clc; clearvars; format long;

%% Setup paths and parameters
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 1;
countY = 10;
sizeX = 1;
sizeY = 1;
gridSize = [countX, countY];
origin = [0, 0, 0];

order = 5;
lambda = 5;
k = 2 * pi / lambda;

%% Create mesh and related objects
myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin, 0, 0, 0);
myMeshCluster = mom.Mesh2DCluster();
myMeshCluster.loadMesh(myRectMesh1.getAsMesh());

quadrature = mom.Quadrature(order);
basisFunctions = mom.BasisFunctions(myMeshCluster);

rooftopBasisFunctions = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize);
rooftopBasisFunctions.generatePoints();

strategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, rooftopBasisFunctions, k);

conjugateGradientSolver = mom.FastImpedanceMatrix.ConjugateGradientSolver(strategyCtx);
conjugateGradientSolver.generateOperators();

feederPosition = ceil(countY / 2);
V = zeros(basisFunctions.getCount(), 1);
V(feederPosition) = 1;
excitationVectorX = [];
excitationVectorY = V;

tolerance = 1e-5;
maxIterationCount = 1000;

[currentX, currentY, converged, residueNorm, iterationCount] = conjugateGradientSolver.solveCurrent(excitationVectorX, excitationVectorY, tolerance, maxIterationCount)
