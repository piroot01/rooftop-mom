clc; clearvars; format long;

% load the source
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

% dipole lengths
dipoleLengths = 100:100:10000;

% times
timeToCompute = zeros(4, numel(dipoleLengths));

countX = 1;
sizeX = 1;
sizeY = 1;

origin = [0, 0, 0];
order = 5;

memoryElementwise = 1e5;
memory_MB = memoryElementwise * 8 / 1e6;

quadrature = mom.Quadrature(order);

for n = 1:numel(dipoleLengths)
    times = zeros(4, 1);
    for m = 1:4
        fprintf("Processing dipole with length %d.\n", dipoleLengths(n));

        lambda = 2 * dipoleLengths(n) * sizeY;
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
        gridSize{1} = [countX, dipoleLengths(n)];

        tic;

        Z = chunkedOptimiser.compute(strategyContext, gridSize);

        time1 = toc;

        times(m) = time1;
    end

    timeToCompute(:, n) = times;
end

disp(timeToCompute);
