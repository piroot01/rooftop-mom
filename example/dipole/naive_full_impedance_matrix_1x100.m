clc; clearvars;

% Include toolbox paths.
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', '..', 'src');
if exist(srcDir, 'dir')
    addpath(srcDir);
else
    warning('Toolbox path not found: %s\nConstructors may fail if the toolbox is not installed there.', srcDir);
end

%% Antenna mesh setup:
% The mesh is described by how many faces are in each direction.
countX = 100;
countY = 1;

% Size of each elementary cell (face) of the mesh.
sizeX = 1;
sizeY = 1;

% Generate a primitive description of the mesh. This is a lightweight
% structure holding only the basic tiling info (no vertices yet).
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

% Convert the primitive into a 3D-oriented rectangle object.
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);

% Generate the actual mesh of the antenna.
mesh = rect.generateMesh();

%% Quadrature / rooftop basis setup:
% Order controls accuracy of Gauss-Legendre quadrature. Higher => more accurate,
% but memory and CPU cost grow quickly (roughly O(N^4) with order N).
order = 5;

% Quadrature for two rectangular faces (Gauss-Legendre).
quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

% Build rooftop basis on the mesh (also computes any auxiliary info needed)
rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quad);
rooftops.generate();

%% Assembly planning and wavenumber:
% How many elements of the distance matrix to compute per chunk.
elementCount = 1e6;

% Planner controls how chunks are arranged during impedance assembly.
planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

% Example wavenumber: lambda = 2 * sizeX * countX (lambda-half dipole),
% so k = 2 * pi / lambda = pi / (sizeX * countX).
k = pi / (sizeX * countX);

% Bundle required data for the assembler.
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

%% Assemble full impedance matrix:
% FullAssembler efficiently composes magnetic and electric parts for us.
fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

% Start assembly (result stored in 'output'; matrix accessible via output.value())
tic;
output = fullAssembler.generateFull(context);
elapsed = toc;
fprintf('Assembly finished in %.2f s.\n', elapsed);
