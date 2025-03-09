clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 2;
countY = 2;
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];
order = 2;

tic;

% 1.step: generate a mesh
myRectMeshPrimitive = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh = mom.MeshRect(myRectMeshPrimitive, origin, 0, 0, 0);
myMesh = myRectMesh.getAsMesh();

% 2. step: create Quadrature
quadrature = mom.Quadrature(order);

% 3. step: create BasisFunctions
basisFuncs = mom.BasisFunctions(myMesh);

% 4. step: generate quadrature points and rooftop points
basisFuncs.generatePoints(myMesh, quadrature);

% 5. step: get the basis function faces
basisFuncs.getFaces()

% qP = basisFunc.getQuadPoints();
% rP = basisFunc.getRooftopPoints();

toc;