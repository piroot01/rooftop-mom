clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 0.4;
sizeY = 0.5;
origin = [0, 0, 0];
order = 2;

tic;

% 1.step: generate a mesh
myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMeshPrimitive2 = mom.MeshRectPrimitive(countX, countY, sizeX * 2, sizeY / 2, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin, 0, 0, 0);
myRectMesh2 = mom.MeshRect(myRectMeshPrimitive2, origin + [0, 0, 2], 0, 0, 0);
myMesh1 = myRectMesh1.getAsMesh();
myMesh2 = myRectMesh2.getAsMesh();

meshCluster = mom.Mesh2DCluster();
meshCluster.loadMesh(myMesh1);
%meshCluster.loadMesh(myMesh2);

% 2. step: create Quadrature
quadrature = mom.Quadrature(order);

% 3. step: create BasisFunctions
basisFunctions = mom.BasisFunctions(meshCluster);

% 4. step: create ComputeStrategyOptimiser
computeStrategyOptimiser = mom.ImpedanceMatrix.ComputeStrategyOptimiser(basisFunctions, [countX, countY]);
[mask, interactions] = computeStrategyOptimiser.getMaskAndInteractions();

disp(mask);
disp(interactions);

% 5. step: create ComputeChunkedContext
memoryElementWise = 48;
memory_MB = memoryElementWise * 8 / 1e6;

basisFunctionCount = basisFunctions.getCount();
computeChunkedContext = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);

tic;

chunks = computeChunkedContext.generateChunksAllInteractions(quadrature.getOrder(), basisFunctionCount);

toc;
tic;

slices = computeChunkedContext.generateSlicesFromChunks(chunks, quadrature.getOrder()^2, basisFunctionCount);

toc;

chunks = computeChunkedContext.generateChunksSpecificInteractions(quadrature.getOrder(), interactions, basisFunctionCount);
celldisp(chunks);

slices = computeChunkedContext.generateSlicesFromChunks(chunks, quadrature.getOrder()^2, basisFunctionCount);
celldisp(slices)
