clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 10;
countY = 10;
gridSize = [countX, countY];
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];

order = 5;
lambda = 5;
k = 2 * pi / lambda;

memoryElementwise = 1e7;
memory_MB = memoryElementwise * 8 / 1e6;

myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin, 0, 0, 0);
myMeshCluster = mom.Mesh2DCluster();
myMeshCluster.loadMesh(myRectMesh1.getAsMesh());

quadrature = mom.Quadrature(order);
basisFunctions = mom.BasisFunctions(myMeshCluster);

rooftops = mom.RooftopBasisFunctions(myMeshCluster, basisFunctions.getFaces());
rooftops.generatePoints(quadrature);

fastRooftops = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize);
fastRooftops.generatePoints();

%
% FFT based impedance matrix solver
%

fastStrategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, fastRooftops, k);

conjugateGradientSolver = mom.FastImpedanceMatrix.ConjugateGradientSolver(fastStrategyCtx);

tic;

conjugateGradientSolver.generateOperators();
Z = conjugateGradientSolver.computeImpedanceMatrix1();

t1 = toc;

fprintf("Number of unknowns = %d\n", length(Z));
fprintf("Time to compute the impedance matrix using FFT solver: %f\n", t1);

%
% optimised regular impedance matrix solver
%

strategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, rooftops, k);
chunkedStrategyCtx = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);

strategy = mom.ImpedanceMatrix.ComputeChunked(chunkedStrategyCtx);

chunkedOptimiser = mom.ImpedanceMatrix.ComputeChunkedOptimiser(strategy, chunkedStrategyCtx);
gridSizeCell{1} = [countX, countY];

tic;

chunkedOptimiser.compute(strategyCtx, gridSizeCell);

t2 = toc;

fprintf("Time to compute the impedance matrix using an optimised regular solver: %f\n", t2);

%
% regular impedance matrix solver
%

tic;

strategy.compute(strategyCtx);

t3 = toc;

fprintf("Time to compute the impedance matrix using a regular solver: %f\n", t3);
