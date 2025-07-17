clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

%
% inputs
%

% quadrature order
order = 5;

% memory per one chunk element-wise
memoryElementwise = 1e7;
memory_MB = memoryElementwise * 8 / 1e6;

% CG settings
tolerance = 1e-4;
maxIterationCount = 20000;

% plane wave parameters
direction = [0, 0, -1];

% polarization parallel to the y axis
polarization = [0, 1, 0];

%
% preparation
%

countX = 3;
countY = 3;
gridSize = [countX, countY];
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];

%
% antenna genome
%

antennaGenome = true(countX, countY);
antennaGenome(2, 2) = false;

disp(antennaGenome);


lambda = 2 * countY * sizeY;
k = 2 * pi / lambda;

% mesh generation
myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin, 0, 0, 0);
myMeshCluster = mom.Mesh2DCluster();
myMeshCluster.loadMesh(myRectMesh1.getAsMesh());

quadrature = mom.Quadrature(order);
basisFunctions = mom.BasisFunctions(myMeshCluster);

fastRooftops = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize);
fastStrategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, fastRooftops, k);
perturbedSolver = mom.FastImpedanceMatrix.PerturbedConjugateGradientSolver(fastStrategyCtx);

%
% operator generation
%

tic;

perturbedSolver.generateOperators(memory_MB, antennaGenome);

timeToComputeOperators = toc;

fprintf("Time to compute the Green's operator: %f\n", timeToComputeOperators);

%
% excitation - plane wave
%

planeWaveExcitatiton = mom.ImpedanceMatrix.ComputePlaneWaveExcitation();

V = planeWaveExcitatiton.compute(fastStrategyCtx, direction, polarization);
V = V(:);

excitationVectorX = V(1:basisFunctions.getCountInDirection('x'));
excitationVectorY = V((basisFunctions.getCountInDirection('x') + 1):end);

%
% PCG-FFT solver
%

tic;

[currentX, currentY, converged, residueNorm, iterationCount] = perturbedSolver.solveCurrentPCG(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);

timeToSolveTheSystem = toc;

fprintf("Time to solve using PCG-FFT: %f\n", timeToSolveTheSystem);
fprintf("Number of iterations: %s\n", mat2str(iterationCount));

% old school

rooftopBasisFunctions = mom.RooftopBasisFunctions(myRectMesh1.getAsMesh(), basisFunctions.getFaces());
rooftopBasisFunctions.generatePoints(quadrature);

strategyContext = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, rooftopBasisFunctions, k);
computeChunkedContext = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);
strategy = mom.ImpedanceMatrix.ComputeChunked(computeChunkedContext);

chunkedOptimiser = mom.ImpedanceMatrix.ComputeChunkedOptimiser(strategy, computeChunkedContext);

clear gridSize;
gridSize{1} = [countX, countY];

tic;

Z = chunkedOptimiser.compute(strategyContext, gridSize);

toc;

I = Z \ V
