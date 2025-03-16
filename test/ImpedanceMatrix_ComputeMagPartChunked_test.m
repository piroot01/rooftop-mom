clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 2;
countY = 2;
sizeX = 1.5;
sizeY = 0.5;
origin = [0, 0, 0];
order = 2;

lambda = 6;
k = 2 * pi / lambda;

memoryElementwise = 100000;
memory_MB = memoryElementwise * 8 / 1e6;

% 1. step: generate a mesh
myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMeshPrimitive2 = mom.MeshRectPrimitive(countX, countY, sizeX * 1, sizeY * 1, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin, 0, 0, 0);
myRectMesh2 = mom.MeshRect(myRectMeshPrimitive2, origin + [0, 0, 1], 0, 0, 0);
myMeshCluster = mom.Mesh2DCluster();
myMeshCluster.loadMesh(myRectMesh1.getAsMesh());
%myMeshCluster.loadMesh(myRectMesh2.getAsMesh());

% 2. step: create Quadrature
quadrature = mom.Quadrature(order);

% 3. step: create BasisFunctions
basisFuncs = mom.BasisFunctions(myMeshCluster);

% 4. step: generate quadrature points and rooftop points
basisFuncs.generatePoints(quadrature);

% 5. step: create a strategy for mag part of Z
stratCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFuncs, k);
magStratContext = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB, order, basisFuncs.getCount());
magStrat = mom.ImpedanceMatrix.ComputeMagPartChunked(magStratContext);
Zm = magStrat.compute(stratCtx);

disp(Zm);