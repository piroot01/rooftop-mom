clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 1;
sizeY = 1;
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

myRectMeshPrimitive1.getPatchCountX();

% 2. step: create Quadrature
quadrature = mom.Quadrature(order);

% 3. step: create BasisFunctions
basisFuncs = mom.BasisFunctions(myMeshCluster);
disp(basisFuncs.getFaces());

% 4. step: generate quadrature points and rooftop points
basisFuncs.generatePoints(quadrature);

compStratOpt = mom.ImpedanceMatrix.ComputeStrategyOptimiser(basisFuncs, [countX, countY]);

tic;

[mask, interactions] = compStratOpt.getMaskAndInteractions();

disp(interactions);

toc;
