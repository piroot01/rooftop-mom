clc; clearvars; format long;

% load the source
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

% The task is to display three graphs:
%   1. The classical algorithm
%   2. The best‐case scenario of the accelerated interaction-symmetries approach
%   3. The worst-case scenario of the same approach
%
% We test dipole counts from 100 to 10000 (i.e., 9901 values).  For the worst case,
% we tile rectangular patches—whose dimensions are chosen by the findRectangle
% function—to best approximate each dipole count, thereby minimizing symmetry.
%
% Each method is executed 4 times; the first run is discarded, and the last three
% runs are used to compute statistics.

function [Mbest, Nbest, fbest, diffbest] = findRectangle(C)
% findRectangle  Find integer M,N such that (M-1)*N + M*(N-1) ≈ C
%   [Mbest, Nbest, fbest, diffbest] = findRectangle(C) returns the dimensions
%   Mbest, Nbest, the resulting sum fbest, and the absolute error
%   diffbest = |fbest - C|.
%
%   Example:
%     [M,N,f,d] = findRectangle(12345);
%     fprintf('Best M=%d, N=%d: f=%d (error %d)\n', M,N,f,d);

    % Input validation
    if C <= 0 || ~isscalar(C) || C~=floor(C)
        error('C must be a single positive integer.');
    end
    
    % Bound on M from N>=2:
    %   (C + M)/(2M-1) >= 2  ⇒  C >= 3M - 2  ⇒  M <= (C+2)/3
    Mmax = max(2, floor((C + 2)/3));
    
    % Initialize best
    diffbest = inf;
    Mbest = 2;
    Nbest = 2;
    fbest = (Mbest-1)*Nbest + Mbest*(Nbest-1);
    
    % Search over M
    for M = 2:Mmax
        realN = (C + M) / (2*M - 1);
        for N = [floor(realN), ceil(realN)]
            if N < 2, continue; end
            f = (M-1)*N + M*(N-1);
            d = abs(f - C);
            if d < diffbest
                diffbest = d;
                Mbest    = M;
                Nbest    = N;
                fbest    = f;
                if d == 0
                    return;  % exact match found
                end
            end
        end
    end
end

% dipole lengths
dipoleLengths = 100:100:10000;

% times
timeToCompute = zeros(4, numel(dipoleLengths));

sizeX = 1;
sizeY = 1;

origin = [0, 0, 0];
order = 3;

memoryElementwise = 1e6;
memory_MB = memoryElementwise * 8 / 1e6;

quadrature = mom.Quadrature(order);

testCount = numel(dipoleLengths);
times = zeros(testCount, 3);
basisFunctionCounts = zeros(testCount, 2);

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

    rooftopBasisFunctions = mom.RooftopBasisFunctions(dipoleMesh, basisFunctions.getFaces());

    rooftopBasisFunctions.generatePoints(quadrature);

    strategyContext = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, rooftopBasisFunctions, k);
    computeChunkedContext = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);
    strategy = mom.ImpedanceMatrix.ComputeChunked(computeChunkedContext);

    chunkedOptimiser = mom.ImpedanceMatrix.ComputeChunkedOptimiser(strategy, computeChunkedContext);
    gridSize{1} = [countX, countY];

    [countXRect, countYRect, ~, ~] = findRectangle(countY);
    lambdaRect = 2 * countYRect * sizeY;
    kRect = 2 * pi / lambdaRect;

    rectPrimitive = mom.MeshRectPrimitive(countXRect, countYRect, sizeX, sizeY, []);
    rectPrimitiveTransformed = mom.MeshRect(rectPrimitive, origin, 0, 0, 0);
    rectMesh = rectPrimitiveTransformed.getAsMesh();

    basisFunctionsRect = mom.BasisFunctions(rectMesh);

    rooftopBasisFunctionsRect = mom.RooftopBasisFunctions(rectMesh, basisFunctionsRect.getFaces());

    rooftopBasisFunctionsRect.generatePoints(quadrature);

    strategyContextRect = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctionsRect, rooftopBasisFunctionsRect, kRect);
    computeChunkedContextRect = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);
    strategyRect = mom.ImpedanceMatrix.ComputeChunked(computeChunkedContextRect);

    chunkedOptimiserRect = mom.ImpedanceMatrix.ComputeChunkedOptimiser(strategyRect, computeChunkedContextRect);
    gridSizeRect{1} = [countXRect, countYRect];

    fprintf("Dipole length: %d, Rect size: %d x %d\n", dipoleLengths(n), countXRect, countYRect);
    fprintf("Number of rooftops: %d, %d\n", basisFunctions.getCount(), basisFunctionsRect.getCount());
    basisFunctionCounts(n, :) = [basisFunctions.getCount(), basisFunctionsRect.getCount()];

    for m = 1:4
        % dipole using the classsical approach
        if (sizeY < 1100)
            tic;

            strategyRect.compute(strategyContext);

            time_classical = toc;
        else
            time_classical = nan;
        end

        fprintf("[%d/%d]: Classical approach finished %f\n", m, 4, time_classical);

        % accelerated approach computation - best case
        tic;

        chunkedOptimiser.compute(strategyContext, gridSize);

        time_accelerated_best = toc;

        fprintf("[%d/%d]: Accelerated approach best case finished %f\n", m, 4, time_accelerated_best);

        % accelerated approach computation - worst case
        tic;

        chunkedOptimiserRect.compute(strategyContextRect, gridSizeRect);

        time_accelerated_worst = toc;

        fprintf("[%d/%d]: Accelerated approach worst case finished %f\n", m, 4, time_accelerated_worst);

        % store the times
        timesPerSize(m, :) = [time_classical, time_accelerated_best, time_accelerated_worst];
    end

    time_classical_avg = sum(timesPerSize(2:4, 1)) / 3;
    time_accelerated_best_avg = sum(timesPerSize(2:4, 2)) / 3;
    time_accelerated_worst_avg = sum(timesPerSize(2:4, 3)) / 3;

    fprintf("Times for a run [%d/%d]: %f, %f, %f\n\n", n, numel(dipoleLengths), time_classical_avg, time_accelerated_best_avg, time_accelerated_worst_avg);
    times(n, :) = [time_classical_avg, time_accelerated_best_avg, time_accelerated_worst_avg];
end

%%

% for tikz
for k = 1:36
    fprintf("            (%f, %f)\n", basisFunctionCounts(k, 1), times(k, 1));
end