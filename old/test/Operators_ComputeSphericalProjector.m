clc; clearvars; close all; format long;

% load the source
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 2;
countY = 2;
sizeX = 1;
sizeY = 1;

origin = [0, 0, 0];

% quadrature order
order = 2;

lambda = 2 * countY * sizeY;
k = 2 * pi / lambda;

% memory given to the mom
memoryElementwise = 10e5;
memory_MB = memoryElementwise * 8 / 1e6;

% 1. step: generate a mesh
dipolePrimitive = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
% rotate the dipole such that the dipole is in the z axis
dipolePrimitiveTransformed = mom.MeshRect(dipolePrimitive, origin, pi/2, 0, 0);
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

direction = [1, 1, 0];
polarization = [1, 1, 0];
planeWaveExcitatiton = mom.ImpedanceMatrix.ComputePlaneWaveExcitation();

V = planeWaveExcitatiton.compute(strategyContext, direction, polarization);

%feederPosition = ceil(basisFunctions.getCount() / 2);
%V = zeros(1, basisFunctions.getCount());
%V(feederPosition) = 1;

I = Z / V;

L = 1;

S = mom.Operators.ComputeSphericalProjector.compute(strategyContext, k, L)
