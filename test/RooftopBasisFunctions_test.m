clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 2;
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
meshCluster.loadMesh(myMesh2);

%meshCluster.getFaces()

% 2. step: create Quadrature
quadrature = mom.Quadrature(order);

% 3. step: create BasisFunctions
basisFunctions = mom.BasisFunctions(meshCluster);

basisFunctions.getFaces()
%basisFunctions.getFacesReindexed()

% 4. step: create RooftopBasisFunctions
%rooftopBasisFunctions = mom.RooftopBasisFunctions(meshCluster, basisFunctions.getFaces());
%rooftopBasisFunctions.generatePoints(quadrature);
