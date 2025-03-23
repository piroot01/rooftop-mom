clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 2;
countY = 2;
sizeX = 1.5;
sizeY = 0.5;
origin = [0, 0, 0];
order = 5;

lambda = 5;
k = 2 * pi / lambda;

memoryElementwise = 10000000;
memory_MB = memoryElementwise * 8 / 1e6;

% 1. step: generate a mesh
myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMeshPrimitive2 = mom.MeshRectPrimitive(countX, countY, sizeX * 1.3, sizeY * 0.3, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin, 0, 0, 0);
myRectMesh2 = mom.MeshRect(myRectMeshPrimitive2, origin + [0, 0, 1], 0, 0, 0);
myMeshCluster = mom.Mesh2DCluster();
myMeshCluster.loadMesh(myRectMesh1.getAsMesh());
myMeshCluster.loadMesh(myRectMesh2.getAsMesh());

% 2. step: create Quadrature
quadrature = mom.Quadrature(order);

% 3. step: create BasisFunctions
basisFuncs = mom.BasisFunctions(myMeshCluster);

% 4. step: create RooftopBasisFunctions
rooftops = mom.RooftopBasisFunctions(myMeshCluster, basisFuncs.getFaces());

% 5. step: generate quadrature points and rooftop points
rooftops.generatePoints(quadrature);

% 5. step: create a strategy for mag part of Z
stratCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFuncs, rooftops, k);
magStratContext = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);
magStrat = mom.ImpedanceMatrix.ComputeMagPartChunked(magStratContext);

tic;

Zm = magStrat.compute(stratCtx);

toc;

disp(Zm);

%disp(size(Zm));
