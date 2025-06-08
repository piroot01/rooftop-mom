clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

% constants
c0 = 299792458; % m/s

%
% inputs
%

kaMin = 10e-2;
kaMax = 1000;

numberOfRooftopsPerOneWavelength = 10;
rectangleSideRatio = 2;

patchCountPerShorterSide = ceil((kaMax * numberOfRooftopsPerOneWavelength) / (pi * sqrt(1 + floor(rectangleSideRatio)^2)));

% the main input of the proram
ka = 1000;

% quadrature order
order = 3;

% memory per one chunk element-wise
memoryElementwise = 1e7;
memory_MB = memoryElementwise * 8 / 1e6;

% CG settings
tolerance = 1e-4;
maxIterationCount = 10000;

%
% preparation
%

% the shorter side will be aligned with the x axis
countX = patchCountPerShorterSide;
countY = floor(rectangleSideRatio) * patchCountPerShorterSide;

fprintf("Size of the radiator: %d x %d\n", countX, countY);

gridSize = [countX, countY];

% the patch dimensios stays 1 m x 1 m
sizeX = 1;
sizeY = 1;

origin = [0, 0, 0];

% calculation of the lambda from the ka
lambda = ((patchCountPerShorterSide * pi) / ka) * sqrt(1 + floor(rectangleSideRatio)^2)
k = 2 * pi / lambda;

% shift to center the mesh (probably not neccessery)
shift = [countX/2, countY/2, 0];

% mesh generation
myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin - shift, 0, 0, 0);
myMeshCluster = mom.Mesh2DCluster();
myMeshCluster.loadMesh(myRectMesh1.getAsMesh());

quadrature = mom.Quadrature(order);
basisFunctions = mom.BasisFunctions(myMeshCluster);

fprintf("Unknown count: %d\n", basisFunctions.getCount());

fastRooftops = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize);

fastStrategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, fastRooftops, k);

conjugateGradientSolver = mom.FastImpedanceMatrix.ConjugateGradientSolver(fastStrategyCtx);

%
% operator generation
%

tic;

conjugateGradientSolver.generateOperators(memory_MB);

timeToComputeOperators = toc;

fprintf("Time to compute the Green's operator: %f\n", timeToComputeOperators);

%
% excitation - plane wave perpendicular to the radiator
%

direction = [0, 0, -1];

% polarization parallel to the y axis
polarization = [0, 1, 0];


planeWaveExcitatiton = mom.ImpedanceMatrix.ComputePlaneWaveExcitation();

V = planeWaveExcitatiton.compute(fastStrategyCtx, direction, polarization);
V = V(:);

excitationVectorX = V(1:basisFunctions.getCountInDirection('x'));
excitationVectorY = V((basisFunctions.getCountInDirection('x') + 1):end);

%
% PCG-FFT solver
%

fprintf("Starting the solver ...\n");

tic;

[currentX, currentY, converged, residueNorm, iterationCount] = conjugateGradientSolver.solveCurrentPCG(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);

timeToSolveTheSystem = toc;

fprintf("Time to solve using PCG-FFT: %f\n", timeToSolveTheSystem);
fprintf("Number of iterations: %s\n", mat2str(iterationCount));

%
% Effective area computation
%

% compute the radiation resistance matrix
Z = conjugateGradientSolver.computeImpedanceMatrix1();
Rr = real(Z);

rHat = [0, 0, -1];

I = [currentX(:); currentY(:)];

U = mom.Operators.ComputeRadiationIntensity.computeInDirection(fastStrategyCtx, I, k, rHat)

Pr = 1/2 * real(I' * Rr * I)

G = 4 * pi * U / Pr

Aeff = G * lambda^2 / (4 * pi)

AeffNorm = Aeff / (countX * countY)
