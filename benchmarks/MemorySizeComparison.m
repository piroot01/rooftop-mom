clc; clearvars; format long;

% load the source
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

% The task is to display three graphs for solving the MoM:
%   1. The best-case of the acceleretad interaction-symmetries approach
%   2. The convolutional approach
%
% We test dipole counts from 100 to 10000 (i.e., 9901 values).
%
% Each method is executed 4 times; the first run is discarded, and the last three
% runs are used to compute statistics.

% dipole lengths
%dipoleLengths = logspace(log10(500000), log10(1000000), 10);
dipoleLengths = 500000;

sizeX = 1;
sizeY = 1;

origin = [0, 0, 0];
order = 3;

memoryElementwise = 1e5;
memory_MB = memoryElementwise * 8 / 1e6;

quadrature = mom.Quadrature(order);

testCount = numel(dipoleLengths);
sizes_MiB = zeros(testCount, 2);
basisFunctionCounts = zeros(testCount, 1);

profile on -memory

for n = 1:testCount
    % same for each local run
    fprintf("Processing %d/%d.\n", n, testCount);

    countX = 1;
    countY = dipoleLengths(n);
    lambda = 2 * countY * sizeY;
    k = 2 * pi / lambda;

    dipolePrimitive = mom.MeshRectPrimitive(countX, dipoleLengths(n), sizeX, sizeY, []);
    dipolePrimitiveTransformed = mom.MeshRect(dipolePrimitive, origin, 0, 0, 0);
    dipoleMesh = dipolePrimitiveTransformed.getAsMesh();

    basisFunctions = mom.BasisFunctions(dipoleMesh);

    basisFunctionCounts(n) = basisFunctions.getCount();

    fprintf("Processing a dipole of length: %d\n", basisFunctionCounts(n));

    %rooftopBasisFunctions = mom.RooftopBasisFunctions(dipoleMesh, basisFunctions.getFaces());

    %rooftopBasisFunctions.generatePoints(quadrature);

    %strategyContext = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, rooftopBasisFunctions, k);
    %computeChunkedContext = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);
    %strategy = mom.ImpedanceMatrix.ComputeChunked(computeChunkedContext);

    %chunkedOptimiser = mom.ImpedanceMatrix.ComputeChunkedOptimiser(strategy, computeChunkedContext);
    gridSize{1} = [countX, countY];

    fastRooftops = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize{1});
    fastRooftops.generatePoints();
    
    fastStrategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, fastRooftops, k);
    
    conjugateGradientSolver = mom.FastImpedanceMatrix.ConjugateGradientSolver(fastStrategyCtx);

    conjugateGradientSolver.generateOperators(memory_MB);

    %Z = conjugateGradientSolver.computeImpedanceMatrix2();

    %Z_MiB = whos('Z').bytes / 2^20;
    Z_MiB = nan;

    [~, operators_MiB] = conjugateGradientSolver.getOperatorSizeMiB();

    sizes_MiB(n, :) = [Z_MiB, operators_MiB];

    fprintf("Sizes for a run [%d/%d]: %f, %f\n\n", n, numel(dipoleLengths), Z_MiB, operators_MiB);
end

profile off

profile viewer

%%

% for tikz
for k = 1:length(basisFunctionCounts)
    fprintf("        (%f, %f)\n", basisFunctionCounts(k), sizes_MiB(k, 2));
end