clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 2;
countY = 2;
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];
order = 2;

memoryElementwise = 1024;
memory_MB = memoryElementwise * 8 / 1e6;

% 1. step: generate a mesh
myRectMeshPrimitive = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh = mom.MeshRect(myRectMeshPrimitive, origin, 0, 0, 0);
myMesh = myRectMesh.getAsMesh();

% 2. step: create Quadrature
quadrature = mom.Quadrature(order);

% 3. step: create BasisFunctions
basisFuncs = mom.BasisFunctions(myMesh);

% 4. step: generate quadrature points and rooftop points
basisFuncs.generatePoints(myMesh, quadrature);
basisFuncs.getFaces()

% 5. step: create ImpedanceMatrixChunked
impedMatrix = mom.ImpedanceMatrixChunked(quadrature, basisFuncs, memory_MB);

% 6. step: compute the impedance matrix
impedMatrix.compute();