clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

% only dipole antennas will make sense
countX = 4;
countY = 2;
gridSize = [countX, countY];
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];

order = 2;
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

%q1 = rooftops.getQuadraturePoints();

fastRooftops = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize);

%q2 = fastRooftops.getQuadraturePoints(false);

%max(max(max(q1 - q2)))

fastStrategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, fastRooftops, k);

conjugateGradientSolver = mom.FastImpedanceMatrix.ConjugateGradientSolver(fastStrategyCtx);

conjugateGradientSolver.generateOperators(memory_MB);

% for very large matrices
%tic;

%Z = conjugateGradientSolver.computeImpedanceMatrix2();

%toc;

strategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, rooftops, k);
chunkedStrategyCtx = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);

strategy = mom.ImpedanceMatrix.ComputeChunked(chunkedStrategyCtx);

chunkedOptimiser = mom.ImpedanceMatrix.ComputeChunkedOptimiser(strategy, chunkedStrategyCtx);
gridSizeCell{1} = [countX, countY];

Z = chunkedOptimiser.compute(strategyCtx, gridSizeCell);

fprintf("Number of unknowns = %d\n", basisFunctions.getCount());

direction = [0, 0, -1];
polarization = [0, 1, 0];
planeWaveExcitatiton = mom.ImpedanceMatrix.ComputePlaneWaveExcitation();

V = planeWaveExcitatiton.compute(fastStrategyCtx, direction, polarization);
V = V(:);

%feederPosition = ceil(basisFunctions.getCount() / 2);
%V = zeros(basisFunctions.getCount(), 1);
%V(feederPosition) = 1;

%
% CG-FFT solver
%

excitationVectorX = [];
excitationVectorY = V;

restart = 100;
tolerance = 1e-10;
maxIterationCount = 10000;

tic;

%[currentX, currentY, converged, residueNorm, iterationCount] = conjugateGradientSolver.solveCurrentPGMRES(excitationVectorX, excitationVectorY, restart, tolerance, maxIterationCount);
[currentX, currentY, converged, residueNorm, iterationCount] = conjugateGradientSolver.solveCurrentPCG(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);
%[currentX, currentY, converged, residueNorm, iterationCount] = conjugateGradientSolver.solveCurrentBiCGSTAB(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);
%[currentX, currentY, converged, residueNorm, iterationCount] = conjugateGradientSolver.solveCurrentCG(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);

t1 = toc;

fprintf("Time to solve using CG-FFT: %f\n", t1);
fprintf("Number of iterations: %s\n", mat2str(iterationCount));

current = [currentX(:); currentY(:)];

%
% build-in matlab solver
%

tic;

I = Z \ V;

t2 = toc;

fprintf("Time to solve using build-in '\\' solver: %f\n", t2);

norm(I)
max(current - I)
