clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

% only dipole antennas will make sense
countX = 300;
countY = 33;
gridSize = [countX, countY];
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];

order = 3;
lambda = 2 * countY * sizeY;
k = 2 * pi / lambda;

memoryElementwise = 1e5;
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

fastStrategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, fastRooftops, k);

conjugateGradientSolver = mom.FastImpedanceMatrix.ConjugateGradientSolver(fastStrategyCtx);

conjugateGradientSolver.generateOperators();

% for very large matrices
%Z = conjugateGradientSolver.computeImpedanceMatrix1();

%{

strategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, rooftops, k);
chunkedStrategyCtx = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);

strategy = mom.ImpedanceMatrix.ComputeChunked(chunkedStrategyCtx);

chunkedOptimiser = mom.ImpedanceMatrix.ComputeChunkedOptimiser(strategy, chunkedStrategyCtx);
gridSizeCell{1} = [countX, countY];

Z = chunkedOptimiser.compute(strategyCtx, gridSizeCell);

%}

fprintf("Number of unknowns = %d\n", basisFunctions.getCount());

feederPosition = ceil(basisFunctions.getCount() / 2);
V = zeros(basisFunctions.getCount(), 1);
V(feederPosition) = 1;

%
% CG-FFT solver
%

excitationVectorX = [];
excitationVectorY = V;

restart = 100;
tolerance = 1e-3;
maxIterationCount = 2000;

tic;

%[currentX, currentY, converged, residueNorm, iterationCount] = conjugateGradientSolver.solveCurrentPGMRES(excitationVectorX, excitationVectorY, restart, tolerance, maxIterationCount);
[currentX, currentY, converged, residueNorm, iterationCount] = conjugateGradientSolver.solveCurrentPCG(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);

t1 = toc;

fprintf("Time to solve using CG-FFT: %f\n", t1);
fprintf("Number of iterations: %s\n", mat2str(iterationCount));

fprintf("currentY norm = %f\n", norm(currentY));

%{

current = [currentX(:); currentY(:)];

Zin = (current' * Z * current) / (abs(current(feederPosition))^2);

fprintf("Input impedance of the dipole Zin = %s Ohms.\n", num2str(Zin));

%
% build-in matlab solver
%

tic;

I = Z \ V;

t2 = toc;

fprintf("Time to solve using build-in '\\' solver: %f\n", t2);

Zin = (I' * Z * I) / (abs(I(feederPosition))^2);

fprintf("Input impedance of the dipole Zin = %s Ohms.\n", num2str(Zin));

%}
