clc; clearvars; format long;

% load the source
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

% dipole dimensions
countX = 1;
countY = 10;
sizeX = 1;
sizeY = 1;

% dipole origin
origin = [0, 0, 0];

% quadrature order
order = 5;

% wavelength and wavenumber - here we choose a lambda = dipole length * 2
%lambda = 2 * countY * sizeY;
lambda = 5;
k = 2 * pi / lambda;

% memory given to the mom
memoryElementwise = 10e6;
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

% now let's compute input impedance for a feeder placed in the middle of the dipole
feederPosition = ceil(basisFunctions.getCount() / 2);
V = zeros(basisFunctions.getCount(), 1);
V(feederPosition) = 1;

tic;

I = Z \ V

toc;

norm(I)

Zin = (I' * Z * I) / (abs(I(feederPosition))^2);

fprintf("Input impedance of the dipole Zin = %s Ohms.\n", num2str(Zin));
