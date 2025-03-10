clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 100;
countY = 1;
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
meshCluster.loadMesh(myMesh2);

% meshCluster.getFaces()

% 2. step: create Quadrature
quadrature = mom.Quadrature(order);

% 3. step: create BasisFunctions
basisFuncs = mom.BasisFunctions(meshCluster);

% 4. step: generate quadrature points and rooftop points
basisFuncs.generatePoints(quadrature);

% 5. step: get the basis function faces
basisFuncs.getFaces();

basisFuncs.getRooftopLengths();
% qP = basisFunc.getQuadPoints();
% rP = basisFunc.getRooftopPoints();

% basisFuncs.getFuncIDs()

toc;