clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 20;
countY = 20;
sizeX = 0.43;
sizeY = 1.2;
origin = [0, 0, 0];
order = 5;
gridSize = [countX, countY];

lambda = 5;
k = 2 * pi / lambda;

memoryElementwise = uint32(10000000);
memory_MB = memoryElementwise * 8 / 1e6;

% 1. step: generate a mesh
myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMeshPrimitive2 = mom.MeshRectPrimitive(countX, countY, sizeX * 1.3, sizeY * 0.3, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin, 0, 0, 0);
myRectMesh2 = mom.MeshRect(myRectMeshPrimitive2, origin + [0, 0, 1], 0, 0, 0);
myMeshCluster = mom.Mesh2DCluster();
myMeshCluster.loadMesh(myRectMesh1.getAsMesh());
%myMeshCluster.loadMesh(myRectMesh2.getAsMesh());

% 2. step: create Quadrature
quadrature = mom.Quadrature(order);

% 3. step: create BasisFunctions
basisFuncs = mom.BasisFunctions(myRectMesh1.getAsMesh());

basisFuncs.getFacesY()

% 4. step: create RooftopBasisFunctions
rooftops = mom.RooftopBasisFunctions(myRectMesh1.getAsMesh(), basisFuncs.getFaces());

% 5. step: generate quadrature points and rooftop points
rooftops.generatePoints(quadrature);

% 5. step: create a strategy for mag part of Z
stratCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFuncs, rooftops, k);
magStratContext = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);

% 6. step: create a strategy and optimiser
magStrat = mom.ImpedanceMatrix.ComputeMagPartChunked(magStratContext);
chunkedOptimiser = mom.ImpedanceMatrix.ComputeChunkedOptimiser(magStrat, magStratContext);

tic;

Zm = magStrat.compute(stratCtx);

toc;

tic;

ZmOptimiser = chunkedOptimiser.compute(stratCtx, gridSize);

toc;

%Zm(11:22, 11)
%ZmOptimiser(11:22, 11)

all(abs(Zm - ZmOptimiser) < 1e7, 'all')

%[~, ~, idx] = unique(round(Zm, 6), 'stable');
%mask1 = reshape(idx, size(Zm))

%[~, ~, idx] = unique(round(ZmOptimiser, 6), 'stable');
%mask2 = reshape(idx, size(ZmOptimiser))

%mask1 - mask2
