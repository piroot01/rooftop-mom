clc; clearvars; format long;

% load the source
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

% The task is to display three graphs for constructing the impedance matrix:
%   1. The best-case of the acceleretad interaction-symmetries approach
%   2. The convolutional approach
%
% We test dipole counts from 100 to 10000 (i.e., 9901 values).
%
% Each method is executed 4 times; the first run is discarded, and the last three
% runs are used to compute statistics.

% dipole lengths
dipoleLengths = logspace(log10(14000), log10(17000), 3);
%dipoleLengths = 20000;

% times
timeToCompute = zeros(4, numel(dipoleLengths));

sizeX = 1;
sizeY = 1;

origin = [0, 0, 0];
order = 3;

memoryElementwise = 1e5;
memory_MB = memoryElementwise * 8 / 1e6;

quadrature = mom.Quadrature(order);

testCount = numel(dipoleLengths);
times = zeros(testCount, 3);
basisFunctionCounts = zeros(testCount, 1);

for n = 1:testCount
    timesPerSize = zeros(4, 3);
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

    rooftopBasisFunctions = mom.RooftopBasisFunctions(dipoleMesh, basisFunctions.getFaces());

    rooftopBasisFunctions.generatePoints(quadrature);

    strategyContext = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, rooftopBasisFunctions, k);
    computeChunkedContext = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);
    strategy = mom.ImpedanceMatrix.ComputeChunked(computeChunkedContext);

    chunkedOptimiser = mom.ImpedanceMatrix.ComputeChunkedOptimiser(strategy, computeChunkedContext);
    gridSize{1} = [countX, countY];

    fastRooftops = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize{1});
    fastRooftops.generatePoints();
    
    fastStrategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, fastRooftops, k);
    
    conjugateGradientSolver = mom.FastImpedanceMatrix.ConjugateGradientSolver(fastStrategyCtx);

    for m = 1:4
        % convolutional approach prepare operators
        %{
        tic;

            conjugateGradientSolver.generateOperators();

        time_conv_prep = toc;

        fprintf("[%d/%d]: Convolutional approach finished preparing operators in %f\n", m, 4, time_conv_prep);

        % convolutional approach construct the impedance matrix
        tic;
            
            Z = conjugateGradientSolver.computeImpedanceMatrix2();
        
        time_conv_imped = toc;

        fprintf("[%d/%d]: Convolutional approach finished in %f\n", m, 4, time_conv_imped);
        %}
        % interaction symmetries approach computation - best case
        tic;

        Z = chunkedOptimiser.compute(strategyContext, gridSize);

        time_accelerated_best = toc;

        clear Z;

        fprintf("[%d/%d]: Accelerated approach best case finished %f\n", m, 4, time_accelerated_best);
        
        % store the times
        timesPerSize(m, :) = [nan, nan, time_accelerated_best];
    end

    time_conv_prep_avg = sum(timesPerSize(2:4, 1)) / (length(timesPerSize) - 1);
    time_conv_imped_avg = sum(timesPerSize(2:4, 2)) / (length(timesPerSize) - 1);
    time_accelerated_best_avg = sum(timesPerSize(2:4, 3)) / (length(timesPerSize) - 1);

    fprintf("Times for a run [%d/%d]: %f, %f, %f\n\n", n, numel(dipoleLengths), time_conv_prep_avg, time_conv_imped_avg, time_accelerated_best_avg);
    times(n, :) = [time_conv_prep_avg, time_conv_imped_avg, time_accelerated_best_avg];
end

%%

% for tikz
for k = 1:length(basisFunctionCounts)
    fprintf("        (%f, %f)\n", basisFunctionCounts(k), times(k, 3));
end