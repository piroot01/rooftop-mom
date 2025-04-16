clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 1;
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];

order = 3;
lambda = 5;
k = 2 * pi / lambda;

sizes = 100:100:20000;
times = zeros(2, numel(sizes));

n = 1;
for countY = sizes
    gridSize = [countX, countY];
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

    feederPosition = ceil(countY / 2);
    V = zeros(basisFunctions.getCount(), 1);
    V(feederPosition) = 1;
    excitationVectorX = [];
    excitationVectorY = V;

    tolerance = 1e-3;
    maxIterationCount = 5000;

    tic;

    [currentX, currentY, converged, residueNorm, iterationCount] = conjugateGradientSolver.solveCurrent(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);

    times(2, n) = toc;

    fprintf("[%d/%d]: currentYNorm = %d, iterationCount = %d\n", n, numel(sizes), norm(currentY), iterationCount);
    n = n + 1;
end
