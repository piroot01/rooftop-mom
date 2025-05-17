clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 1;
countY = 1000;
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];
order = 3;

lambda = 5;
k = 2 * pi / lambda;

memoryElementwise = 10e5;
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

% 4. step: create RooftopBasisFunctions
rooftops = mom.RooftopBasisFunctions(myMeshCluster, basisFuncs.getFaces());

% 5. step: generate quadrature points and rooftop points
rooftops.generatePoints(quadrature);

stratCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFuncs, rooftops, k);
compStratCtx = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);

%tic;

%magStrat = mom.ImpedanceMatrix.ComputeMagPartChunked(compStratCtx);
%elStrat = mom.ImpedanceMatrix.ComputeElPartChunked(compStratCtx);
%Zm = magStrat.compute(stratCtx);
%Ze = elStrat.compute(stratCtx);

%toc;

%disp(Zm - Ze);

tic;

combStrat = mom.ImpedanceMatrix.ComputeChunked(compStratCtx);
combStrat.compute(stratCtx);

toc;
