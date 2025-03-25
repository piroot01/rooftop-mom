clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 100;
countY = 1;
sizeX = 1;
sizeY = 0.6;
origin = [0, 0, 0];
order = 2;

gridSizes{1} = [countX, countY];
gridSizes{2} = [countX - 1, countY];
gridSizes{3} = [countX - 1, countY];

lambda = 5;
k = 2 * pi / lambda;

memoryElementwise = 10000;
memory_MB = memoryElementwise * 8 / 1e6;

% 1. step: generate a mesh
myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMeshPrimitive2 = mom.MeshRectPrimitive(countX - 1, countY, sizeX * 2, sizeY * 2, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin, 0, 0, 0);
myRectMesh2 = mom.MeshRect(myRectMeshPrimitive2, origin + [0, 0, 1], 0, 0, 0);
myRectMesh3 = mom.MeshRect(myRectMeshPrimitive2, origin + [0, 0, -2], 0, 0, 0);
myMeshCluster = mom.Mesh2DCluster();
myMeshCluster.loadMesh(myRectMesh1.getAsMesh());
myMeshCluster.loadMesh(myRectMesh2.getAsMesh());
myMeshCluster.loadMesh(myRectMesh3.getAsMesh());

% 2. step: create Quadrature
quadrature = mom.Quadrature(order);

% 3. step: create BasisFunctions
basisFuncs = mom.BasisFunctions(myMeshCluster);

% 4. step: create RooftopBasisFunctions
rooftops = mom.RooftopBasisFunctions(myMeshCluster, basisFuncs.getFaces());

% 5. step: generate quadrature points and rooftop points
rooftops.generatePoints(quadrature);

% 6. step: create a strategy for mag part of Z
stratCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFuncs, rooftops, k);
magStratContext = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);

% 7. step: create a strategy and optimiser
magStrat = mom.ImpedanceMatrix.ComputeMagPartChunked(magStratContext);
chunkedOptimiser = mom.ImpedanceMatrix.ComputeChunkedOptimiser(magStrat, magStratContext);

%tic;

%Zm = magStrat.compute(stratCtx);

%toc;

tic;

ZmOptimiser = chunkedOptimiser.compute(stratCtx, gridSizes);

toc;

disp(size(ZmOptimiser));

%Zm(11:22, 11)
%ZmOptimiser(11:22, 11)

%all(abs(Zm - ZmOptimiser) < 1e7, 'all')

%[~, ~, idx] = unique(round(Zm, 6), 'stable');
%mask1 = reshape(idx, size(Zm))

%[~, ~, idx] = unique(round(ZmOptimiser, 6), 'stable');
%mask2 = reshape(idx, size(ZmOptimiser))

%mask1 - mask2
