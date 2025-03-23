clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 1;
countY = 3;
sizeX = 0.4;
sizeY = 0.5;
origin = [0, 0, 0];
order = 2;

tic;

% 1.step: generate a mesh
myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMeshPrimitive2 = mom.MeshRectPrimitive(countX + 1, countY, sizeX * 2, sizeY / 2, []);
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
% mesh
%basisFuncs = mom.BasisFunctions(myMesh1);
%basisFuncs.getFaces()
%basisFuncs.getCount()
%celldisp(basisFuncs.getFacesAsCell());
%celldisp(basisFuncs.getFacesXAsCell());
%basisFuncs.getIDs()

% cluster
basisFuncsCluster = mom.BasisFunctions(meshCluster);
celldisp(basisFuncsCluster.getFacesAsCell())
basisFuncsCluster.getFaces()
basisFuncsCluster.getGridSizes()
%basisFuncsCluster.getCount()
%celldisp(basisFuncsCluster.getFacesAsCell());
%celldisp(basisFuncsCluster.getFacesXAsCell());
%basisFuncsCluster.getIDs()

toc;
