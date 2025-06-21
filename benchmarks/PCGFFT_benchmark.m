clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

%profile clear;
%profile -memory on;

%
% inputs
%

% unknown number
N = 2000;

% quadrature order
order = 5;

% memory per one chunk element-wise
memoryElementwise = 1e7;
memory_MB = memoryElementwise * 8 / 1e6;

% CG settings
tolerance = 1e-4;
maxIterationCount = 20000;

%
% preparation
%

function [countX, countY] = generateSizeForUnknownCount(N)
    countX = 1;
    countY = N + 1;
end

[countX, countY] = generateSizeForUnknownCount(N);
gridSize = [countX, countY];
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];

% l = lambda / 2
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

conjugateGradientSolver = mom.FastImpedanceMatrix.ConjugateGradientSolver(fastStrategyCtx);

%
% operator generation
%

tic;

conjugateGradientSolver.generateOperators(memory_MB);

timeToComputeOperators = toc;

fprintf("Time to compute the Green's operator: %f\n", timeToComputeOperators);

%
% excitation - middle delta-gap feeder
%

%feederPosition = ceil(basisFunctions.getCount() / 2);
%V = zeros(basisFunctions.getCount(), 1);
%V(feederPosition) = 1;

direction = [0, 0, -1];

% polarization parallel to the y axis
polarization = [0, 1, 0];

% rooftops suitable for the plane wave excitation
% TODO: unite with the fastRooftops
rooftopBasisFunctions = mom.RooftopBasisFunctions(myRectMesh1.getAsMesh(), basisFunctions.getFaces());
rooftopBasisFunctions.generatePoints(quadrature);

excitationStrategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, rooftopBasisFunctions, k);

planeWaveExcitatiton = mom.ImpedanceMatrix.ComputePlaneWaveExcitation();

V = planeWaveExcitatiton.compute(excitationStrategyCtx, direction, polarization);

excitationVectorX = V(1:basisFunctions.getCountInDirection('x'));
excitationVectorY = V((basisFunctions.getCountInDirection('x') + 1):end);

%
% PCG-FFT solver
%

tic;

[currentX, currentY, converged, residueNorm, iterationCount] = conjugateGradientSolver.solveCurrentPCG(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);

timeToSolveTheSystem = toc;

fprintf("Time to solve using CG-FFT: %f\n", timeToSolveTheSystem);
fprintf("Number of iterations: %s\n", mat2str(iterationCount));

%profile off;
%profile viewer;



