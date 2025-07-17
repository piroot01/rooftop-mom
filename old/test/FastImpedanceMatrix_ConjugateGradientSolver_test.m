clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 1;
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];

order = 5;

sizes = 100:100:5000;
times = zeros(3, numel(sizes));

n = 1;
for countY = sizes
    gridSize = [countX, countY];

    lambda = 2 * countY * sizeY;
    k = 2 * pi / lambda;

    myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
    myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin, 0, 0, 0);
    myMeshCluster = mom.Mesh2DCluster();
    myMeshCluster.loadMesh(myRectMesh1.getAsMesh());

    quadrature = mom.Quadrature(order);
    basisFunctions = mom.BasisFunctions(myMeshCluster);

    rooftopBasisFunctions = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize);
    rooftopBasisFunctions.generatePoints();

    strategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, rooftopBasisFunctions, k);

    conjugateGradientSolver = mom.FastImpedanceMatrix.ConjugateGradientSolver(strategyCtx);

    tic;

    conjugateGradientSolver.generateOperators();

    times(1, n) = toc;

    feederPosition = ceil(basisFunctions.getCount() / 2);
    V = zeros(basisFunctions.getCount(), 1);
    V(feederPosition) = 1;

    excitationVectorX = [];
    excitationVectorY = V;

    tolerance = 1e-4;
    maxIterationCount = 5000;

    timesPerIteration = zeros(2, 4);

    for k = 1:4
        tic;

        [~, ~, ~, ~, iterationCount] = conjugateGradientSolver.solveCurrent(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);

        timesPerIteration(1, k) = toc;

        tic;

        [~, ~, ~, ~, iterationCountPrecond] = conjugateGradientSolver.solveCurrentPrecond(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);

        timesPerIteration(2, k) = toc;
    end

    times(2, n) = sum(timesPerIteration(1, 2:end)) / 3;
    times(3, n) = sum(timesPerIteration(2, 2:end)) / 3;

    fprintf("[%d/%d] [%d]: iterationCount = %d, iterationCountPrecond = %d\n", n, numel(sizes), length(V),  iterationCount, iterationCountPrecond);
    n = n + 1;
end
