clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 2;
countY = 2;
sizeX = 0.5;
sizeY = 1.5;
origin = [0, 0, 0];
order = 2;

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

disp(myRectMesh1.getAsMesh().vertices);

% 2. step: create Quadrature
quadrature = mom.Quadrature(order);

% 3. step: create BasisFunctions
basisFuncs = mom.BasisFunctions(myMeshCluster);
disp(basisFuncs.getFaces());

% 4. step: generate quadrature points and rooftop points
basisFuncs.generatePoints(quadrature);

% 5. step: create a strategy for mag part of Z
stratCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFuncs);
elStratContext = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB, order, basisFuncs.getCount());
elStrat = mom.ImpedanceMatrix.ComputeElPartChunked(elStratContext);
Ze = elStrat.compute(stratCtx);

disp(Ze);