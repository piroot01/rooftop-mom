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
%dipoleLengths = logspace(log10(12000), log10(300000), 5);
dipoleLengths = 10000000;

% iteration counts for the iterative solvers
iterationCounts = zeros(numel(dipoleLengths), 2);

sizeX = 1;
sizeY = 1;

origin = [0, 0, 0];
order = 3;

memoryElementwise = 1e7;
memory_MB = memoryElementwise * 8 / 1e6;

quadrature = mom.Quadrature(order);

testCount = numel(dipoleLengths);
times = zeros(testCount, 4);
basisFunctionCounts = zeros(testCount, 1);


for n = 1:testCount
    timesPerSize = zeros(4, 4);
    iterationCountsPerSize = zeros(4, 2);
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

    feederPosition = ceil(basisFunctions.getCount() / 2);
    V = zeros(basisFunctions.getCount(), 1);
    V(feederPosition) = 1;

    excitationVectorX = [];
    excitationVectorY = V;

    tolerance = 1e-5;
    maxIterationCount = 10000;

    %Z = chunkedOptimiser.compute(strategyContext, gridSize);

    for m = 1:4
        % convolutional approach prepare operators
        tic;

            conjugateGradientSolver.generateOperators(memory_MB);

        time_conv_prep = toc;

        fprintf("[%d/%d]: Convolutional approach finished preparing operators in %f\n", m, 4, time_conv_prep);

        % convolutional approach - CG
        %{
        tic;
            
            [~, ~, ~, ~, iterationCountCG] = conjugateGradientSolver.solveCurrentCG(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);
        
        time_conv_cg = toc;

        fprintf("[%d/%d]: Convolutional - CG - approach finished in %f, iterations: %d\n", m, 4, time_conv_cg, iterationCountCG);
        %}

        % convolutional approach - PCG
        tic;
            
            [~, ~, ~, ~, iterationCountPCG] = conjugateGradientSolver.solveCurrentPCG(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);
        
        time_conv_pcg = toc;

        fprintf("[%d/%d]: Convolutional - PCG - approach finished in %f, iterations: %d\n", m, 4, time_conv_pcg, iterationCountPCG);

        % build-in matlab solver
        %{
        tic;

        I = Z \ V;

        time_accelerated_solve = toc;

        fprintf("[%d/%d]: Accelerated approach best case finished %f\n", m, 4, time_accelerated_solve);
        %}
        
        % store the times
        timesPerSize(m, :) = [time_conv_prep, nan, time_conv_pcg, nan];
        iterationCountsPerSize(m, :) = [nan, iterationCountPCG];
    end

    time_conv_prep_avg = sum(timesPerSize(2:4, 1)) / (length(timesPerSize) - 1);
    time_conv_cg_avg = sum(timesPerSize(2:4, 2)) / (length(timesPerSize) - 1);
    time_conv_pcg_avg = sum(timesPerSize(2:4, 3)) / (length(timesPerSize) - 1);
    time_accelerated_solve_avg = sum(timesPerSize(2:4, 4)) / (length(timesPerSize) - 1);

    iterationsCount_cg_avg = sum(iterationCountsPerSize(2:4, 1)) / 3;
    iterationsCount_pcg_avg = sum(iterationCountsPerSize(2:4, 2)) / 3;

    fprintf("Times for a run [%d/%d]: %f, %f, %f, %f\n\n", n, numel(dipoleLengths), time_conv_prep_avg, time_conv_cg_avg, time_conv_pcg_avg, time_accelerated_solve_avg);
    times(n, :) = [time_conv_prep_avg, time_conv_cg_avg, time_conv_pcg_avg, time_accelerated_solve_avg];
    iterationCounts(n, :) = [iterationsCount_cg_avg, iterationsCount_pcg_avg];
end

%%

% for tikz
for k = 1:length(basisFunctionCounts)
    fprintf("        (%f, %f)\n", basisFunctionCounts(k), times(k, 1));
end