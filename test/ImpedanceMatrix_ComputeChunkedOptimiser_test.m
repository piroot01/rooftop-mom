clc; clearvars; format long;

% load the source
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 10;
countY = 11;
sizeX = 1;
sizeY = 1;

origin = [0, 0, 0];

% quadrature order
order = 5;

lambda = 5;
k = 2 * pi / lambda;

% memory given to the mom
memoryElementwise = 10e5;
memory_MB = memoryElementwise * 8 / 1e6;

% 1. step: generate a mesh
dipolePrimitive = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
dipolePrimitiveTransformed = mom.MeshRect(dipolePrimitive, origin, 0, 0, 0);
dipoleMesh = dipolePrimitiveTransformed.getAsMesh();

% 2. step: create Quadrature
quadrature = mom.Quadrature(order);

% 3. step: create BasisFunctions
basisFunctions = mom.BasisFunctions(dipoleMesh);

% 4. step: create RooftopBasisFunctions
rooftopBasisFunctions = mom.RooftopBasisFunctions(dipoleMesh, basisFunctions.getFaces());

% 5. step: generate quadrature points and rooftop points
rooftopBasisFunctions.generatePoints(quadrature);

% 6. step: crate a supporting objects for the Z computation
strategyContext = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, rooftopBasisFunctions, k);
computeChunkedContext = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);
strategy = mom.ImpedanceMatrix.ComputeChunked(computeChunkedContext);

% 7. step: initialize the optimiser
chunkedOptimiser = mom.ImpedanceMatrix.ComputeChunkedOptimiser(strategy, computeChunkedContext);
gridSize{1} = [countX, countY];

tic;

Z = chunkedOptimiser.compute(strategyContext, gridSize);

toc;

disp(size(Z));
